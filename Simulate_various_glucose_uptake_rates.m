%% Simulations for a range of glucose uptake rates

% Timing: ~  s

% No saturation factor is used here.

% Simulated results will be saved in the folder 'Results'.

tic;
load('pcLactis_Model.mat');
model = pcLactis_Model;

%% Optimization setting.

rxnID = 'R_dummy_assumed_Monomer';
osenseStr = 'Maximize';

%% Parameters.
GAM = 42;%ATP coefficient in the new biomass equation.
NGAM = 2.5; %(mmol/gCDW/h)
f_unmodeled = 0.42; %proportion of unmodeled protein in total protein (g/g)

model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,f] = ChangeUnmodeledProtein(model,f_unmodeled);

kcat_glc = 180;%kcat value of glucose transporter

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
model = changeRxnBounds(model,'R_M_PROTS_LLA',0,'b');
model = changeRxnBounds(model,'R_M_PROTS_LLA_v2',0,'b');
model = changeRxnBounds(model,'R_M_PROTS_LLA_v3',0,'b');
model = changeRxnBounds(model,'R_M_MGt2pp_rvs',0,'b');%block infinite h[e]

%% Search max mu.

glc_range = -2:-10:-32;
fluxes_simulated = zeros(length(model.rxns),length(glc_range));

factor_k = 1;%global saturation factor
f_transporter = 1;%no constraint on glucose transporter
factor_glc = 1;%assume glucose transporter is fully saturated

for i = 1:length(glc_range)
    glc_ex = glc_range(i);
    model = changeRxnBounds(model,'R_M_EX_glc_LPAREN_e_RPAREN_',glc_ex,'l');
    
    mu_low = 0;
    mu_high = 0.8;
    
    while mu_high-mu_low > 0.01
        mu_mid = (mu_low+mu_high)/2;
        disp(['q_glc = ' num2str(glc_ex) '; mu = ' num2str(mu_mid)]);
        model = changeRxnBounds(model,'R_biomass_dilution',mu_mid,'b');
        model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*mu_mid,'l');

        fileName = WriteLPSatFactor(model,mu_mid,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,factor_glc,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,solME_status,~] = ReadSoplexResult(fileName_out,model);
        
        if strcmp(solME_status,'optimal')
            mu_low = mu_mid;
        else
            mu_high = mu_mid;
        end
    end
    
    model = changeRxnBounds(model,'R_biomass_dilution',mu_low,'b');
    model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*mu_low,'l');
	
	fileName = WriteLPSatFactor(model,mu_mid,f,osenseStr,rxnID,factor_k,...
                                f_transporter,kcat_glc,factor_glc,...
                                Info_enzyme,...
                                Info_mRNA,...
                                Info_protein,...
                                Info_ribosome,...
                                Info_tRNA);
	command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
	system(command,'-echo');
    fileName_out = 'Simulation.lp.out';
    [~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model);
    
    if strcmp(solME_status,'optimal')
        fluxes_simulated(:,i) = solME_full;
    else
        fluxes_simulated(:,i) = zeros(length(model.rxns),1);
    end
end

cd Results/;
save('Svgur_fluxes.mat','fluxes_simulated');
cd ../;

clear;
toc;