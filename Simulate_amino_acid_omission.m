%% Simulate amino acid omission
% Timing: ~ 8000 s

% Simulated results will be saved in the folder 'Results'.

tic;
load('pcLactis_Model.mat');
model = pcLactis_Model;

%% Optimization setting.

rxnID = 'R_dummy_assumed_Monomer';
osenseStr = 'Maximize';

%% Parameters.
GAM = 36;%ATP coefficient in the new biomass equation.
NGAM = 3; %(mmol/gCDW/h)
f_unmodeled = 0.4; %proportion of unmodeled protein in total protein (g/g)

model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,f] = ChangeUnmodeledProtein(model,f_unmodeled);

kcat_glc = 180;%kcat value of glucose transporter
f_transporter = 0.0083;%fraction of glucose transporter in total proteome
factor_k = 1;

%% Data import.
load('Info_enzyme.mat');
load('Info_mRNA.mat');
load('Info_protein.mat');
load('Info_ribosome.mat');
load('Info_tRNA.mat');
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

factor_glc = 1; % unlimited glucose

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

fluxes = zeros(length(model.rxns),length(aa_list)+1);
mu_list = zeros(1,length(aa_list)+1);
    
% Determine amino acid uptake and growth rate in reference
model_ref = model;
mu_low = 0;
mu_high = 0.8;
while mu_high-mu_low > 0.001
    mu_mid = (mu_low+mu_high)/2;
    disp(['Ref: mu = ' num2str(mu_mid)]);
    model_ref = changeRxnBounds(model_ref,'R_biomass_dilution',mu_mid,'b');
    model_ref = changeRxnBounds(model_ref,Exchange_AAs,LBfactor_AAs*mu_mid,'l');
    fileName = WriteLPSatFactor(model_ref,mu_mid,f,osenseStr,rxnID,factor_k,...
        f_transporter,kcat_glc,factor_glc,...
        Info_enzyme,...
        Info_mRNA,...
        Info_protein,...
        Info_ribosome,...
        Info_tRNA);
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
mu_list(1,1) = flux_tmp(strcmp(model_ref.rxns,'R_biomass_dilution'),1);
fluxes(:,1) = flux_tmp;

for j = 1:length(aa_list)
    model_tmp = model;
    aaid = aa_list(j);
    aarxnid = Exchange_AAs(contains(Exchange_AAs,aaid));
    
    mu_low = 0;
    mu_high = 0.8;
    flux_tmp = zeros(length(model.rxns),1);
    while mu_high-mu_low > 0.001
        mu_mid = (mu_low+mu_high)/2;
        disp([cell2mat(aaid) ': mu = ' num2str(mu_mid)]);
        model_tmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_mid,'b');
        model_tmp = changeRxnBounds(model_tmp,Exchange_AAs,LBfactor_AAs*mu_mid,'l');
        model_tmp = changeRxnBounds(model_tmp,aarxnid,0,'l');
        fileName = WriteLPSatFactor(model_tmp,mu_mid,f,osenseStr,rxnID,factor_k,...
            f_transporter,kcat_glc,factor_glc,...
            Info_enzyme,...
            Info_mRNA,...
            Info_protein,...
            Info_ribosome,...
            Info_tRNA);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-12 -o1e-12 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model_tmp);
        if strcmp(solME_status,'optimal')
            mu_low = mu_mid;
            flux_tmp = solME_full;
        else
            mu_high = mu_mid;
        end
    end
    mu_list(1,1+j) = flux_tmp(strcmp(model_tmp.rxns,'R_biomass_dilution'),1);
    fluxes(:,1+j) = flux_tmp;
end

cd Results/;
save('Saao_fluxes.mat','fluxes');
save('Saao_result.mat','mu_list');
cd ../;

clear;
toc;
