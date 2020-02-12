%% Estimate a global saturation factor 2.

% Timing: ~ 6200 s

% Consider glucose transporter

% Most of the enzymes in the model show increased saturation factor with
% increased growth rate, so we estimate a global saturation factor for
% them. For the rest, we assume that their saturation factor are unchanged.

% Simulated results will be saved in the folder 'Results'.

tic;
load('pcLactis_Model.mat');
model = pcLactis_Model;

%% Optimization setting.
rxnID = 'R_dummy_assumed_Monomer';
% rxnID = 'R_M_3MBALDt_rvs_Enzyme';
osenseStr = 'Maximize';

%% Parameters.
GAM = 36; %ATP coefficient in the new biomass equation.
NGAM = 2; %(mmol/gCDW/h)
f_unmodeled = 0.4; %proportion of unmodeled protein in total protein (g/g)

model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,f] = ChangeUnmodeledProtein(model,f_unmodeled);

kcat_glc = 180;%kcat value of glucose transporter
f_transporter = 0.009;%fraction of glucose transporter in total proteome

%% Data import.
load('Info_enzyme.mat');
load('Info_mRNA.mat');
load('Info_protein.mat');
load('Info_ribosome.mat');
load('Info_tRNA.mat');

%% Set reaction rates.
% Block uptake direction of all the exchange reactions.
idx = cellfun(@(x) contains(x,'R_M_EX_'),model.rxns,'UniformOutput',false);
M_exchange_reactions = model.rxns(cell2mat(idx));
model = changeRxnBounds(model,M_exchange_reactions,0,'l');
clear idx M_exchange_reactions;

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

%% Loop for dilution rate of 0.15 0.3 0.45 0.5 and 0.6.
[~, ~, exchange_raw] = xlsread('Exchange_reaction_setting.xlsx','Exp_bounds');

Exchange_reactions = exchange_raw(2:end,1);
header = exchange_raw(1,1:21);

mu_list = [0.15 0.3 0.45 0.5 0.6];

global_saturation_factor_list = [];

fluxes_global_saturation_factor_changed = [];

for i = 1:length(mu_list)
    mu = mu_list(i);
    
    model_tmp = changeRxnBounds(model,'R_biomass_dilution',mu,'b');

    % Set bounds for some exchange reactions.
    idx = find(contains(header,num2str(mu)));
    replicate = length(idx)/2;
    
    for j = 1:replicate
        lb_idx = idx(2*j-1);
        ub_idx = idx(2*j);
        LB = cell2mat(exchange_raw(2:end,lb_idx));
        UB = cell2mat(exchange_raw(2:end,ub_idx));
        model_tmp = changeRxnBounds(model_tmp,Exchange_reactions,LB,'l');
        model_tmp = changeRxnBounds(model_tmp,Exchange_reactions,UB,'u');

        disp(['mu = ' num2str(mu) '; replicate = ' num2str(j)]);
        
        [factor_k_new,...
         solME_full_new] = SearchSaturationFactor(model_tmp,mu,f,osenseStr,rxnID,...
                                                  f_transporter,kcat_glc,...
                                                  Info_enzyme,...
                                                  Info_mRNA,...
                                                  Info_protein,...
                                                  Info_ribosome,...
                                                  Info_tRNA,...
                                                  2);
                                                
        global_saturation_factor_list = [global_saturation_factor_list;mu factor_k_new];
        
        fluxes_global_saturation_factor_changed = [fluxes_global_saturation_factor_changed solME_full_new];
    end
end

cd Results/;
save('Egsf2_result.mat','global_saturation_factor_list');
save('Egsf2_fluxes.mat','fluxes_global_saturation_factor_changed');
cd ../;
clear;

toc;
