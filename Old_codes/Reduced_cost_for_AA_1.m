%% Reduced cost for AA uptake
% Timing: ~ 53000 s

% Use absolute increase, e.g., 0.1mu mmol/gCDW/h;
% Use unchanged uptake rates for other AAs;
% Use changed saturation factor, catalytic rates, degradation constants...

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
f_transporter = 0.009;%fraction of glucose transporter in total proteome

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
% load glucose concentration of reference state
load('Sglc2_result_with_sf.mat');
selected_points = [1;5;9;13];
glc_conc_list = glc_conc_with_sf(selected_points,:);
clear glc_conc_with_sf;

AA_factor = 0.1; % increase = mu * AA_factor (mmol/gCDW/h)

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

% obtain the global saturation factor
load('Egsf2_result.mat');
x = global_saturation_factor_list(:,1);
y = global_saturation_factor_list(:,2);
% x = x(y ~= 1);
% y = y(y ~= 1);
x = x(~isnan(y));
y = y(~isnan(y));
sf_coeff = x\y;
clear x y;

fluxes_rcAA = zeros(length(model.rxns),length(selected_points)*length(aa_list));
result_rcAA.row = {'mu_ref';'mu';'aa_ref';'aa_bound';'aa_real'};
result_rcAA.column = cell(1,length(selected_points)*length(aa_list));
result_rcAA.data = zeros(length(result_rcAA.row),length(selected_points)*length(aa_list));

for i = 1:length(selected_points)
    
    mu_ref_setting = glc_conc_list(i,1);
	glc_conc = glc_conc_list(i,2);
    factor_glc = glc_conc / (glc_conc + Km);
    increase = mu_ref_setting * AA_factor;
    
    % Determine amino acid uptake in reference
    model_ref = model;
	mu_low = 0;
	mu_high = 0.8;
	while mu_high-mu_low > 0.00001
        mu_mid = (mu_low+mu_high)/2;
        disp(['Ref: Glucose concentration = ' num2str(glc_conc) ' uM; mu = ' num2str(mu_mid)]);
        model_ref = changeRxnBounds(model_ref,'R_biomass_dilution',mu_mid,'b');
        model_ref = changeRxnBounds(model_ref,Exchange_AAs,LBfactor_AAs*mu_mid,'l');
        factor_k = sf_coeff * mu_mid;
        if factor_k > 1
            factor_k = 1;
        end
        fileName = WriteLPSatFactor(model_ref,mu_mid,f,osenseStr,rxnID,factor_k,...
                                        f_transporter,kcat_glc,factor_glc,...
                                        Info_enzyme,...
                                        Info_mRNA,...
                                        Info_protein,...
                                        Info_ribosome,...
                                        Info_tRNA);
        command = sprintf('/Users/cheyu/build/bin/soplex -t1000 -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,solME_status,~] = ReadSoplexResult(fileName_out,model_ref);
        if strcmp(solME_status,'optimal')
            mu_low = mu_mid;
        else
            mu_high = mu_mid;
        end
	end
    
	model_ref = changeRxnBounds(model_ref,'R_biomass_dilution',mu_low,'b');
	model_ref = changeRxnBounds(model_ref,Exchange_AAs,LBfactor_AAs*mu_low,'l');
	factor_k = sf_coeff * mu_low;
	if factor_k > 1
        factor_k = 1;
	end
	fileName = WriteLPSatFactor(model_ref,mu_low,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,factor_glc,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);
	command = sprintf('/Users/cheyu/build/bin/soplex -t1000 -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
	system(command,'-echo');
	fileName_out = 'Simulation.lp.out';
	[~,~,solME_full] = ReadSoplexResult(fileName_out,model_ref);
    
	mu_ref = solME_full(strcmp(model_ref.rxns,'R_biomass_dilution'),1);
    [~, idx_tmp] = ismember(Exchange_AAs,model_ref.rxns);
    q_AAs_ref = solME_full(idx_tmp);
    q_AAs_ref(q_AAs_ref > 0) = 0;
    
    for j = 1:length(aa_list)
        model_tmp = model;
        model_tmp = changeRxnBounds(model_tmp,Exchange_AAs,q_AAs_ref,'l');
        
        aaid = aa_list(j);
        aarxnid = Exchange_AAs(contains(Exchange_AAs,aaid));
        aaref = q_AAs_ref(contains(Exchange_AAs,aaid));
        aalb = aaref - increase;
        model_tmp = changeRxnBounds(model_tmp,aarxnid,aalb,'l');

        mu_low = 0;
        mu_high = 0.8;
    
        while mu_high-mu_low > 0.00001
            mu_mid = (mu_low+mu_high)/2;
            disp(['Glucose concentration = ' num2str(glc_conc) '; ' cell2mat(aaid) '; mu = ' num2str(mu_mid)]);
            model_tmptmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_mid,'b');
            
            factor_k = sf_coeff * mu_mid;
            if factor_k > 1
                factor_k = 1;
            end
            fileName = WriteLPSatFactor(model_tmptmp,mu_mid,f,osenseStr,rxnID,factor_k,...
                                        f_transporter,kcat_glc,factor_glc,...
                                        Info_enzyme,...
                                        Info_mRNA,...
                                        Info_protein,...
                                        Info_ribosome,...
                                        Info_tRNA);
            command = sprintf('/Users/cheyu/build/bin/soplex -t1000 -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
            system(command,'-echo');
            fileName_out = 'Simulation.lp.out';
            [~,solME_status,~] = ReadSoplexResult(fileName_out,model_tmptmp);
            if strcmp(solME_status,'optimal')
                mu_low = mu_mid;
            else
                mu_high = mu_mid;
            end
        end
        
        model_tmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_low,'b');
        factor_k = sf_coeff * mu_low;
        if factor_k > 1
            factor_k = 1;
        end
        fileName = WriteLPSatFactor(model_tmp,mu_low,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,factor_glc,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);
        command = sprintf('/Users/cheyu/build/bin/soplex -t1000 -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model_tmp);
        
        result_rcAA.column(1,(i-1)*20+j) = strcat(num2str(mu_ref_setting),'_',aaid);
        if strcmp(solME_status,'optimal')
            fluxes_rcAA(:,(i-1)*20+j) = solME_full;
            res_tmp = [mu_ref;mu_low;-aaref;-aalb;-solME_full(strcmp(model_tmp.rxns,aarxnid),1)];
            result_rcAA.data(:,(i-1)*20+j) = res_tmp;
        else
            fluxes_rcAA(:,(i-1)*20+j) = zeros(length(model_tmp.rxns),1);
            res_tmp = [mu_ref;0;-aaref;-aalb;0];
            result_rcAA.data(:,(i-1)*20+j) = res_tmp;
        end
    end
end

cd Results/;
save('RcAA_fluxes.mat','fluxes_rcAA');
save('RcAA_result.mat','result_rcAA');
cd ../;

clear;
toc;

%% Figures
load('RcAA_result.mat');

increase_mu = round(result_rcAA.data(2,:),5)-round(result_rcAA.data(1,:),5);
increase_mu(increase_mu < 0) = 0;
increase_aa = round(result_rcAA.data(5,:),5)-round(result_rcAA.data(3,:),5);
increase_aa(increase_aa < 0) = 0;
reduced_cost = increase_mu./increase_aa;
reduced_cost = round(reduced_cost,5);

aaidlist = result_rcAA.column(1:20);
aaidlist = cellfun(@(x) x(5:end),aaidlist,'UniformOutput',false);
c = categorical(aaidlist);

figure('Name','1');
subplot(4,1,1);
bar(c,reduced_cost(1:20));
ylim([0 0.15]);
set(gca,'FontSize',12,'FontName','Helvetica');
title('Original mu = 0.1 /h','FontSize',14,'FontName','Helvetica');
ylabel('Reduced cost','FontSize',14,'FontName','Helvetica');

subplot(4,1,2);
bar(c,reduced_cost(21:40));
ylim([0 0.15]);
set(gca,'FontSize',12,'FontName','Helvetica');
title('Original mu = 0.3 /h','FontSize',14,'FontName','Helvetica');
ylabel('Reduced cost','FontSize',14,'FontName','Helvetica');

subplot(4,1,3);
bar(c,reduced_cost(41:60));
ylim([0 0.15]);
set(gca,'FontSize',12,'FontName','Helvetica');
title('Original mu = 0.5 /h','FontSize',14,'FontName','Helvetica');
ylabel('Reduced cost','FontSize',14,'FontName','Helvetica');

subplot(4,1,4);
bar(c,reduced_cost(61:80));
ylim([0 0.15]);
set(gca,'FontSize',12,'FontName','Helvetica');
title('Original mu = 0.7 /h','FontSize',14,'FontName','Helvetica');
ylabel('Reduced cost','FontSize',14,'FontName','Helvetica');

clear;
