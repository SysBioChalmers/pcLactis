%% FVA for glucose uptake rate 
% Calculate UB and LB of the glucose uptake rate for each glucose
% concentration.

% Timing: ~ 1700 s

% Simulated results will be saved in the folder 'Results'.

% Figures can be obtained by running the codes starting from line 149.

tic;
load('pcLactis_Model.mat');
model = pcLactis_Model;

%% Parameters.
rxnID = 'R_M_EX_glc__D_e';
GAM = 36;%ATP coefficient in the new biomass equation.
NGAM = 3; %(mmol/gCDW/h)
f_unmodeled = 0.4; %proportion of unmodeled protein in total protein (g/g)

model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,f] = ChangeUnmodeledProtein(model,f_unmodeled);

kcat_glc = 180;%kcat value of glucose transporter
Km = 21;%Km of glucose transporter, unit: uM (PMID: 30630406)
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
load('Sglc_result.mat');

D_list = glc_conc_without_sf(:,1);%unit: /h
glc_conc_list = glc_conc_without_sf(:,2)*1.000001;%unit: uM

fluxes_simulated = zeros(length(model.rxns),2*length(D_list));

for i = 1:length(D_list)
    
    D = D_list(i);
	model = changeRxnBounds(model,'R_biomass_dilution',D,'b');
	model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*D,'l');
    
	glc_conc = glc_conc_list(i);
    factor_glc = glc_conc / (glc_conc + Km);
    
    % Minimize glucose uptake rate
    osenseStr = 'Maximize';
    disp(['D = ' num2str(D) '; Minimizing glucose uptake']);
    
	fileName = WriteLPSatFactor(model,D,f,osenseStr,rxnID,factor_k,...
                                f_transporter,kcat_glc,factor_glc,...
                                Info_enzyme,...
                                Info_mRNA,...
                                Info_protein,...
                                Info_ribosome,...
                                Info_tRNA);
                   
	command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
	system(command,'-echo');
	fileName_out = 'Simulation.lp.out';
	[~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model);
    if strcmp(solME_status,'optimal')
        fluxes_simulated(:,2*i-1) = solME_full;
    else
        fluxes_simulated(:,2*i-1) = zeros(length(model.rxns),1);
    end
    
    % Maximize glucose uptake rate 
    osenseStr = 'Minimize';
    disp(['D = ' num2str(D) '; Maximizing glucose uptake']);
    
	fileName = WriteLPSatFactor(model,D,f,osenseStr,rxnID,factor_k,...
                                f_transporter,kcat_glc,factor_glc,...
                                Info_enzyme,...
                                Info_mRNA,...
                                Info_protein,...
                                Info_ribosome,...
                                Info_tRNA);
                   
	command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
	system(command,'-echo');
	fileName_out = 'Simulation.lp.out';
	[~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model);
    
    if strcmp(solME_status,'optimal')
        fluxes_simulated(:,2*i) = solME_full;
    else
        fluxes_simulated(:,2*i) = zeros(length(model.rxns),1);
    end
end

cd Results/;
save('Eguv_fluxes.mat','fluxes_simulated');
cd ../;

clear;

toc;


