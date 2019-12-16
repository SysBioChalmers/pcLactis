%% Sensitivity analysis for total modeled proteome.

% Timing: ~  s

% When changing the total modeled proteome, other growth rate-dependent 
% parameters are not changed. Especially, the mass of glucose transporter
% should be constant.

% With the saturation saturation factor are performed.

% Simulated results will be saved in the folder 'Results'.

% Figures can be obtained by running the codes starting from line .

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
f_transporter_ref = 0.009;

model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');

kcat_glc = 180;%kcat value of glucose transporter
Km = 21;%Km of glucose transporter, unit: uM (PMID: 30630406)

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

glc_list = [2 4 6 8 10 20 40 60 80 100 200 1000 10000 100000 1000000];%unit: /uM
modeled_protein_list = (1-f_unmodeled_ref)*(1:0.01:1.05);
res_mp = zeros(length(modeled_protein_list),4,length(glc_list));
fluxes_mp = zeros(length(model.rxns),length(modeled_protein_list)*length(glc_list));

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

for i = 1:length(glc_list)
    glc_conc = glc_list(i);
    factor_glc = glc_conc / (glc_conc + Km);
    
    for j = 1 % determine reference growth rate
        
        model_tmp = model;
        
        f_unmodeled = 1-modeled_protein_list(j);
        [model_tmp,f] = ChangeUnmodeledProtein(model_tmp,f_unmodeled);
        f_transporter = f_transporter_ref;
        
        mu_low = 0;
        mu_high = 1;
        
        while mu_high-mu_low > 0.000000001
            mu_mid = (mu_low+mu_high)/2;
            disp(['Glucose concentration = ' num2str(glc_conc) ' uM; modeled protein = ' num2str((1-f_unmodeled)) '; mu = ' num2str(mu_mid)]);
            model_tmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_mid,'b');
            model_tmp = changeRxnBounds(model_tmp,Exchange_AAs,LBfactor_AAs*mu_mid,'l');
            factor_k = sf_coeff * mu_mid;
            if factor_k > 1
                factor_k = 1;
            end
            fileName = WriteLPSatFactor(model_tmp,mu_mid,f,osenseStr,rxnID,factor_k,...
                                        f_transporter,kcat_glc,factor_glc,...
                                        Info_enzyme,...
                                        Info_mRNA,...
                                        Info_protein,...
                                        Info_ribosome,...
                                        Info_tRNA);
            command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-15 -o1e-15 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
            system(command,'-echo');
            fileName_out = 'Simulation.lp.out';
            [~,solME_status,~] = ReadSoplexResult(fileName_out,model_tmp);
            if strcmp(solME_status,'optimal')
                mu_low = mu_mid;
            else
                mu_high = mu_mid;
            end
        end
        
        model_tmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_low,'b');
        model_tmp = changeRxnBounds(model_tmp,Exchange_AAs,LBfactor_AAs*mu_low,'l');
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
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-15 -o1e-15 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,~,solME_full] = ReadSoplexResult(fileName_out,model_tmp);
        
        mu = solME_full(strcmp(model_tmp.rxns,'R_biomass_dilution'),1);
        glc = -solME_full(strcmp(model_tmp.rxns,'R_M_EX_glc__D_e'),1);
        arg = -solME_full(strcmp(model_tmp.rxns,'R_M_EX_arg__L_e'),1);

        f_modeled_tmp = 1-f_unmodeled;
        
        res_mp(j,:,i) = [f_modeled_tmp mu mu/(glc*180/1000) arg*174/1000/mu];
                                          %g_CDW/g_glucose  g_arginine/g_CDW
        fluxes_mp(:,(i-1)*length(glc_list)+j) = solME_full;
        mu_ref = mu;
    end
    
	for j = 2:length(modeled_protein_list)
        
        model_tmp = model;
        
        f_unmodeled = 1-modeled_protein_list(j);
        [model_tmp,f] = ChangeUnmodeledProtein(model_tmp,f_unmodeled);
        f_transporter = f_transporter_ref;
        
        mu_low = 0;
        mu_high = 1;
        
        while mu_high-mu_low > 0.000000001
            mu_mid = (mu_low+mu_high)/2;
            disp(['Glucose concentration = ' num2str(glc_conc) ' uM; modeled protein = ' num2str((1-f_unmodeled)) '; mu = ' num2str(mu_mid)]);
            model_tmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_mid,'b');
            model_tmp = changeRxnBounds(model_tmp,Exchange_AAs,LBfactor_AAs*mu_ref,'l');
            factor_k = sf_coeff * mu_ref;
            if factor_k > 1
                factor_k = 1;
            end
            fileName = WriteLPSatFactorTmp(model_tmp,mu_mid,f,osenseStr,rxnID,factor_k,...
                                           f_transporter,kcat_glc,factor_glc,...
                                           Info_enzyme,...
                                           Info_mRNA,...
                                           Info_protein,...
                                           Info_ribosome,...
                                           Info_tRNA,...
                                           mu_ref);
            command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-15 -o1e-15 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
            system(command,'-echo');
            fileName_out = 'Simulation.lp.out';
            [~,solME_status,~] = ReadSoplexResult(fileName_out,model_tmp);
            if strcmp(solME_status,'optimal')
                mu_low = mu_mid;
            else
                mu_high = mu_mid;
            end
        end
        
        model_tmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_low,'b');
        model_tmp = changeRxnBounds(model_tmp,Exchange_AAs,LBfactor_AAs*mu_ref,'l');
        factor_k = sf_coeff * mu_ref;
        if factor_k > 1
            factor_k = 1;
        end
        fileName = WriteLPSatFactorTmp(model_tmp,mu_low,f,osenseStr,rxnID,factor_k,...
                                       f_transporter,kcat_glc,factor_glc,...
                                       Info_enzyme,...
                                       Info_mRNA,...
                                       Info_protein,...
                                       Info_ribosome,...
                                       Info_tRNA,...
                                       mu_ref);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-15 -o1e-15 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,~,solME_full] = ReadSoplexResult(fileName_out,model_tmp);
        
        mu = solME_full(strcmp(model_tmp.rxns,'R_biomass_dilution'),1);
        glc = -solME_full(strcmp(model_tmp.rxns,'R_M_EX_glc__D_e'),1);
        arg = -solME_full(strcmp(model_tmp.rxns,'R_M_EX_arg__L_e'),1);

        f_modeled_tmp = 1-f_unmodeled;
        
        res_mp(j,:,i) = [f_modeled_tmp mu mu/(glc*180/1000) arg*174/1000/mu];
                                          %g_CDW/g_glucose  g_arginine/g_CDW
        fluxes_mp(:,(i-1)*length(glc_list)+j) = solME_full;
	end
end

cd Results/;
save('Samp_result.mat','res_mp');
save('Samp_fluxes.mat','fluxes_mp');
cd ../;

clear;
toc;
