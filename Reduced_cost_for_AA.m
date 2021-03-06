%% Reduced cost for AA uptake
% Timing: ~ 115000 s

% Use relative increase, i.e., 1% of the reference value;
% Keep unchanged uptake rates for other AAs as reference;
% Keep unchanged saturation factor, catalytic rates, degradation constants...as reference;

% Simulated results will be saved in the folder 'Results'.

tic;
load('pcLactis_Model.mat');
model = pcLactis_Model;

%% Optimization setting.

rxnID = 'R_dummy_assumed_Monomer';
osenseStr = 'Maximize';

%% Parameters.
GAM = 36;%ATP coefficient in the new biomass equation.
NGAM = 2; %(mmol/gCDW/h)
f_unmodeled = 0.4; %proportion of unmodeled protein in total protein (g/g)

model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,f] = ChangeUnmodeledProtein(model,f_unmodeled);

kcat_glc = 180;%kcat value of glucose transporter
Km = 21;%Km of glucose transporter, unit: uM (PMID: 30630406)
f_transporter = 0.0082;%fraction of glucose transporter in total proteome
factor_k = 1;

%% Data import.
load('Info_enzyme.mat');
load('Info_mRNA.mat');
load('Info_protein.mat');
load('Info_ribosome.mat');
[~, ~, exchange_raw] = xlsread('Exchange_reaction_setting.xlsx','MaxMu');
[~, ~, aa_raw] = xlsread('Exchange_reaction_setting.xlsx','AA_factors');

Exchange_AAs = aa_raw(2:end,1);
LBfactor_AAs = cell2mat(aa_raw(2:end,2));
clear aa_raw;

%% Set reaction rates.
% Block uptake direction of all the exchange reactions.
idx = cellfun(@(x) contains(x,'R_M_EX_'),model.rxns,'UniformOutput',false);
M_exchange_reactions = model.rxns(cell2mat(idx));
model = changeRxnBounds(model,M_exchange_reactions,0,'l');
clear idx M_exchange_reactions;
% Set bounds for some exchange reactions.
Exchange_reactions = exchange_raw(2:end,1);
LB = cell2mat(exchange_raw(2:end,2));
UB = cell2mat(exchange_raw(2:end,3));
model = changeRxnBounds(model,Exchange_reactions,LB,'l');
model = changeRxnBounds(model,Exchange_reactions,UB,'u');
clear exchange_raw Exchange_reactions LB UB;
% Block some reactions in the M model.
model = changeRxnBounds(model,'R_M_biomass_LLA',0,'b');
model = changeRxnBounds(model,'R_M_biomass_LLA_atpm',0,'b');
model = changeRxnBounds(model,'R_M_PROTS_LLA_v3',0,'b');

% Block other glucose transporters
model = changeRxnBounds(model,'R_M_GLCpts_1',0,'b');
model = changeRxnBounds(model,'R_M_GLCt2_fwd',0,'b');

% Block one of ADH isozymes llmg_0955
model = changeRxnBounds(model,'R_M_ALCD2x_1_rvs',0,'b');

% Block pyruvate oxidase
model = changeRxnBounds(model,'R_M_PYROX_1',0,'b');

%% Main part.
% load glucose concentration of reference state
load('Sglc_result.mat');
% selected_points = [1,9,17,21,26];
selected_points = 1:8:25;
glc_conc_list = glc_conc_without_sf(selected_points,:);
clear glc_conc_without_sf;

AA_factor = 0.01; % increase = reference value * AA_factor 

% load AA id
aa_list = {'ala'
           'arg'
           'asn'
           'asp'
           'cys'
           'gln'
           'glu'
           'gly'
           'his'
           'ile'
           'leu'
           'lys'
           'met'
           'phe'
           'pro'
           'ser'
           'thr'
           'trp'
           'tyr'
           'val'};

fluxes_rcAA = zeros(length(model.rxns),length(selected_points)*length(aa_list));
fluxes_ref = zeros(length(model.rxns),length(selected_points));
result_rcAA.row = {'mu_ref';'mu';'aa_ref';'aa_bound';'aa_real'};
result_rcAA.column = cell(1,length(selected_points)*length(aa_list));
result_rcAA.data = zeros(length(result_rcAA.row),length(selected_points)*length(aa_list));

for i = 1:length(selected_points)
    
    mu_ref_setting = glc_conc_list(i,1);
	glc_conc = glc_conc_list(i,2);
    factor_glc = glc_conc / (glc_conc + Km);
    increase = mu_ref_setting * AA_factor;
    
    % Determine amino acid uptake and growth rate in reference
    model_ref = model;
	mu_low = 0;
	mu_high = 0.8;
    while mu_high-mu_low > 0.000000001
        mu_mid = (mu_low+mu_high)/2;
        disp(['Ref: Glucose concentration = ' num2str(glc_conc) ' uM; mu = ' num2str(mu_mid)]);
        model_ref = changeRxnBounds(model_ref,'R_biomass_dilution',mu_mid,'b');
        model_ref = changeRxnBounds(model_ref,Exchange_AAs,LBfactor_AAs*mu_mid,'l');
        fileName = WriteLPSatFactor(model_ref,mu_mid,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,factor_glc,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-12 -o1e-12 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model_ref);
        if strcmp(solME_status,'optimal')
            mu_low = mu_mid;
            flux_tmp = solME_full;
        else
            mu_high = mu_mid;
        end
    end
    
	mu_ref = flux_tmp(strcmp(model_ref.rxns,'R_biomass_dilution'),1);
    [~, idx_tmp] = ismember(Exchange_AAs,model_ref.rxns);
    q_AAs_ref = flux_tmp(idx_tmp);
    q_AAs_ref(q_AAs_ref > 0) = 0;
    q_AAs_ref = q_AAs_ref*1.00001;
    
    % Estimate reference state using the new AA uptakes and growth rate
    model_tmp = model;
    model_tmp = changeRxnBounds(model_tmp,Exchange_AAs,q_AAs_ref,'l');
    disp(['Calculating Ref(new) for Glucose concentration = ' num2str(glc_conc) ' and mu = ' num2str(mu_ref)]);
    model_tmptmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_ref,'b');
    fileName = WriteLPSatFactorTmp(model_tmptmp,mu_ref,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,factor_glc,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    mu_ref);
    command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-12 -o1e-12 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
    system(command,'-echo');
    fileName_out = 'Simulation.lp.out';
    [~,~,solME_full] = ReadSoplexResult(fileName_out,model_tmptmp);
    flux_tmp = solME_full;
    mu_ref_new = flux_tmp(strcmp(model_tmp.rxns,'R_biomass_dilution'),1);
    solME_full_ref = flux_tmp;
    fluxes_ref(:,i) = solME_full_ref;
    [~, idx_tmp] = ismember(Exchange_AAs,model_ref.rxns);
    q_AAs_ref_new = solME_full_ref(idx_tmp);
    q_AAs_ref_new(q_AAs_ref_new > 0) = 0;
    % Reduced cost analysis for each AA
    for j = 1:length(aa_list)
        model_tmp = model;
        model_tmp = changeRxnBounds(model_tmp,Exchange_AAs,q_AAs_ref,'l');
        
        aaid = aa_list(j);
        aarxnid = Exchange_AAs(contains(Exchange_AAs,aaid));
        aaref = q_AAs_ref_new(contains(Exchange_AAs,aaid));
        aalb = aaref * (1 + AA_factor);
        model_tmp = changeRxnBounds(model_tmp,aarxnid,aalb,'l');
        
        mu_low = 0;
        mu_high = 0.8;
        
        while mu_high-mu_low > 0.000000001
            mu_mid = (mu_low+mu_high)/2;
            disp(['Glucose concentration = ' num2str(glc_conc) '; ' cell2mat(aaid) '; mu = ' num2str(mu_mid)]);
            model_tmptmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_mid,'b');
            
            fileName = WriteLPSatFactorTmp(model_tmptmp,mu_mid,f,osenseStr,rxnID,factor_k,...
                                            f_transporter,kcat_glc,factor_glc,...
                                            Info_enzyme,...
                                            Info_mRNA,...
                                            Info_protein,...
                                            Info_ribosome,...
                                            mu_ref_new);
            command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-12 -o1e-12 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
            system(command,'-echo');
            fileName_out = 'Simulation.lp.out';
            [~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model_tmptmp);
            if strcmp(solME_status,'optimal')
                mu_low = mu_mid;
                flux_tmp = solME_full;
            else
                mu_high = mu_mid;
            end
        end
        result_rcAA.column(1,(i-1)*20+j) = strcat(num2str(mu_ref_setting),'_',aaid);
        fluxes_rcAA(:,(i-1)*20+j) = flux_tmp;
        res_tmp = [mu_ref_new;flux_tmp(strcmp(model_tmp.rxns,'R_biomass_dilution'),1);-aaref;-aalb;-flux_tmp(strcmp(model_tmp.rxns,aarxnid),1)];
        result_rcAA.data(:,(i-1)*20+j) = res_tmp;
    end
end

cd Results/;
save('RcAA_fluxes.mat','fluxes_rcAA');
save('RcAA_fluxes_ref.mat','fluxes_ref');
save('RcAA_result.mat','result_rcAA');
cd ../;

clear;
toc;
