%% Simulate glucose-limited chemostats (objective: minimizing glucose concentration)

% Timing: ~ 10000 s

% With the saturation saturation factor is performed.

% Simulated results will be saved in the folder 'Results'.

% Perform two types of simulations:
% 1. without constraining any AA uptakes
% 2. constrain arg uptake rate based on experimental data but set free other
% AA uptakes.

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

% with saturation factor
% obtain the global saturation factor
load('Egsf2_result.mat');
x = global_saturation_factor_list(:,1);
y = global_saturation_factor_list(:,2);
x = x(y ~= 1);
y = y(y ~= 1);
x = x(~isnan(y));
y = y(~isnan(y));
sf_coeff = x\y;
clear x y;

% obtain extracellular glucose concentration
load('Sglc2_result_with_sf.mat');
glc_conc_list = glc_conc_with_sf(:,2);
% glc_conc_list = glc_conc_list([1;3;5;7;9;11;13]);

% 1. without constraining any AA uptakes
fluxes_simulated_without_arg = zeros(length(model.rxns),length(glc_conc_list));

model = changeRxnBounds(model,'R_M_SERD_L',0,'b');
model = changeRxnBounds(model,'R_M_CYSDS',0,'b');

for i = 1:length(glc_conc_list)
    
    glc_conc = glc_conc_list(i);
    factor_glc = glc_conc / (glc_conc + Km);
    
    LBfactor_AAs_tmp = ones(length(LBfactor_AAs),1)*-1000;
    
    model_ref = model;
    model_ref = changeRxnBounds(model_ref,Exchange_AAs,LBfactor_AAs_tmp,'l');
    
	mu_low = 0;
	mu_high = 1;
    
	while mu_high-mu_low > 0.00001
        mu_mid = (mu_low+mu_high)/2;
        disp(['Without arg: Glucose concentration = ' num2str(glc_conc) ' uM; mu = ' num2str(mu_mid)]);
        model_ref = changeRxnBounds(model_ref,'R_biomass_dilution',mu_mid,'b');
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
	[~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model_ref);
    
    if strcmp(solME_status,'optimal')
        fluxes_simulated_without_arg(:,i) = solME_full;
    else
        fluxes_simulated_without_arg(:,i) = zeros(length(model.rxns),1);
    end
end

% 2. constrain arg uptake rate based on experimental data but set free other
% AA uptakes.
fluxes_simulated_with_arg = zeros(length(model.rxns),length(glc_conc_list));
aa_idx = contains(Exchange_AAs,{'ser';'cys'});

for i = 1:length(glc_conc_list)
    
    glc_conc = glc_conc_list(i);
    factor_glc = glc_conc / (glc_conc + Km);
    
    LBfactor_AAs_tmp = ones(length(LBfactor_AAs),1)*-1000;
    
    model_ref = model;
    model_ref = changeRxnBounds(model_ref,Exchange_AAs,LBfactor_AAs_tmp,'l');
    
	mu_low = 0;
	mu_high = 1;
    
	while mu_high-mu_low > 0.00001
        mu_mid = (mu_low+mu_high)/2;
        disp(['With_arg: Glucose concentration = ' num2str(glc_conc) ' uM; mu = ' num2str(mu_mid)]);
        model_ref = changeRxnBounds(model_ref,'R_biomass_dilution',mu_mid,'b');
        model_ref = changeRxnBounds(model_ref,Exchange_AAs(aa_idx),LBfactor_AAs(aa_idx)*mu_mid,'l');
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
    model_ref = changeRxnBounds(model_ref,Exchange_AAs(aa_idx),LBfactor_AAs(aa_idx)*mu_low,'l');
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
	[~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model_ref);
    
    if strcmp(solME_status,'optimal')
        fluxes_simulated_with_arg(:,i) = solME_full;
    else
        fluxes_simulated_with_arg(:,i) = zeros(length(model.rxns),1);
    end
end

cd Results/;
save('Sglc_fluxes_without_arg.mat','fluxes_simulated_without_arg');
save('Sglc_fluxes_with_arg.mat','fluxes_simulated_with_arg');
cd ../;

clear;
toc;

