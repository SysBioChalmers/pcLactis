%% Sensitivity analysis for glucose transporter.

% Timing: ~ 39000 s

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

f_unmodeled_ref = 0.4;
f_transporter_ref = 0.0082;

model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');

kcat_glc = 180;%kcat value of glucose transporter
Km = 21;%Km of glucose transporter, unit: uM (PMID: 30630406)

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
precision = 1e-9;
org_min_mu = 0;
org_max_mu = 1;
factor_k = 1;

% load glucose concentration of reference state
load('Sglc_result.mat');
% selected_points = [1:2:23,25,26];
selected_points = 1:2:25;
glc_list = glc_conc_without_sf(selected_points,2);
clear glc_conc_without_sf selected_points;

% glc_list = [2 4 6 8 10 20 40 60 80 100 200 1000 10000 100000 1000000];%unit: uM
% glucose_transporter_list = f_transporter_ref*(1:0.01:1.05);
% glc_list = [2 5 10 20 50 100 1000 10000];%unit: uM
glucose_transporter_list = f_transporter_ref*[1,1.01];
res_gt = zeros(length(glucose_transporter_list),length(glc_list));
fluxes_gt = zeros(length(model.rxns),length(glucose_transporter_list)*length(glc_list));

for i = 1:length(glc_list)
    glc_conc = glc_list(i);
    factor_glc = glc_conc / (glc_conc + Km);
    
	for k = 1 % determine reference growth rate
        
        model_tmp = model;
        
        f_transporter = glucose_transporter_list(k);
        f_unmodeled = f_unmodeled_ref;
        [model_tmp,f] = ChangeUnmodeledProtein(model_tmp,f_unmodeled);
        
        mu_low = org_min_mu;
        mu_high = org_max_mu;
        
        while mu_high-mu_low > precision
            mu_mid = (mu_low+mu_high)/2;
            disp(['Glucose concentration = ' num2str(glc_conc) '; glucose transporter = ' num2str(f_transporter) '; mu = ' num2str(mu_mid)]);
            model_tmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_mid,'b');
            model_tmp = changeRxnBounds(model_tmp,Exchange_AAs,LBfactor_AAs*mu_mid,'l');
            fileName = WriteLPSatFactor(model_tmp,mu_mid,f,osenseStr,rxnID,factor_k,...
                                        f_transporter,kcat_glc,factor_glc,...
                                        Info_enzyme,...
                                        Info_mRNA,...
                                        Info_protein,...
                                        Info_ribosome);
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
        mu = flux_tmp(strcmp(model_tmp.rxns,'R_biomass_dilution'),1);
        res_gt(k,i) = mu;
        fluxes_gt(:,(i-1)*length(glucose_transporter_list)+k) = flux_tmp;
        mu_ref = mu;
	end
    
	for k = 2:length(glucose_transporter_list)
        
        model_tmp = model;
        
        f_transporter = glucose_transporter_list(k);
        f_unmodeled = f_unmodeled_ref;
        [model_tmp,f] = ChangeUnmodeledProtein(model_tmp,f_unmodeled);
        
        mu_low = org_min_mu;
        mu_high = org_max_mu;
        
        while mu_high-mu_low > precision
            mu_mid = (mu_low+mu_high)/2;
            disp(['Glucose concentration = ' num2str(glc_conc) '; glucose transporter = ' num2str(f_transporter) '; mu = ' num2str(mu_mid)]);
            model_tmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_mid,'b');
            model_tmp = changeRxnBounds(model_tmp,Exchange_AAs,LBfactor_AAs*mu_ref,'l');
            fileName = WriteLPSatFactorTmp(model_tmp,mu_mid,f,osenseStr,rxnID,factor_k,...
                                           f_transporter,kcat_glc,factor_glc,...
                                           Info_enzyme,...
                                           Info_mRNA,...
                                           Info_protein,...
                                           Info_ribosome,...
                                           mu_ref);
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
        mu = flux_tmp(strcmp(model_tmp.rxns,'R_biomass_dilution'),1);
        res_gt(k,i) = mu;
        fluxes_gt(:,(i-1)*length(glucose_transporter_list)+k) = flux_tmp;
	end
end

cd Results/;
save('Sagt_result.mat','res_gt');
save('Sagt_fluxes.mat','fluxes_gt');
cd ../;

toc;
