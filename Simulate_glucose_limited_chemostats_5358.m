%% Simulate glucose-limited chemostats (D=0.53 to 0.58)
% 
% Timing: ~ 10000 s

% Simulated results will be saved in the folder 'Results'.

tic;
load('pcLactis_Model.mat');
model = pcLactis_Model;

%% Optimization setting.

rxnID = 'R_dummy_assumed_Monomer';
osenseStr = 'Maximize';

%% Parameters.
GAM = 36; %ATP coefficient in the new biomass equation.
NGAM = 2; %(mmol/gCDW/h)
f_unmodeled = 0.4; %proportion of unmodeled protein in total protein (g/g)

model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,f] = ChangeUnmodeledProtein(model,f_unmodeled);

kcat_glc = 180;%kcat value of glucose transporter
Km = 21;%Km of glucose transporter, unit: uM (PMID: 30630406)
f_transporter = 0.0082;%fraction of glucose transporter in total proteome

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

D_list = 0.53:0.005:0.58;%unit: /h
% without saturation factor
factor_k = 1;%global saturation factor
fluxes_simulated_without_sf = zeros(length(model.rxns),length(D_list));
glc_conc_without_sf = zeros(length(D_list),2);

for i = 1:length(D_list)
    
    D = D_list(i);
    
    factor_glc_low = 0;
    factor_glc_high = 1;
    
    while factor_glc_high-factor_glc_low > 0.000001
        factor_glc_mid = (factor_glc_low+factor_glc_high)/2;
        disp(['Without sf: D = ' num2str(D) '; factor_glc = ' num2str(factor_glc_mid)]);
        model = changeRxnBounds(model,'R_biomass_dilution',D,'b');
        model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*D,'l');
        
        fileName = WriteLPSatFactor(model,D,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,factor_glc_mid,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome);

        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-17 -o1e-17 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model);
        
        if strcmp(solME_status,'optimal')
            factor_glc_high = factor_glc_mid;
            flux_tmp = solME_full;
        else
            factor_glc_low = factor_glc_mid;
        end
    end
    fluxes_simulated_without_sf(:,i) = flux_tmp;
    glc_conc = Km * factor_glc_high / (1 - factor_glc_high);
    glc_conc_without_sf(i,1) = D;
    glc_conc_without_sf(i,2) = glc_conc;
end

cd Results/;
save('Sglc5358_fluxes.mat','fluxes_simulated_without_sf');
save('Sglc5358_result.mat','glc_conc_without_sf');
cd ../;

clear;
toc;
