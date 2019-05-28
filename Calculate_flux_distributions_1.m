%% Calculate flux distributions using experimental data (PMID: 25828364)

% Timing: ~ 600 s

% Minimizing glucose uptake rate.

% Experimentally measured exchange reaction rates are used as constraints,
% and simulations are performed by minimizing glucose uptake rate without
% and with saturation factor.

% Simulated results will be saved in the folder 'Results'.

tic;
load('pcLactis_Model.mat');
model = pcLactis_Model;

%% Optimization setting.
rxnID = 'R_M_EX_glc_LPAREN_e_RPAREN_';
osenseStr = 'Maximize';

%% Parameters.
GAM = 42;%ATP coefficient in the new biomass equation.
NGAM = 2.5; %(mmol/gCDW/h)
f_unmodeled = 0.42; %proportion of unmodeled protein in total protein (g/g)

model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,f] = ChangeUnmodeledProtein(model,f_unmodeled);

kcat_glc = 180;%kcat value of glucose transporter
f_transporter = 0.01;%fraction of glucose transporter in total proteome

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

% Set NGAM.
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');

% Block some reactions in the M model.
model = changeRxnBounds(model,'R_M_biomass_LLA',0,'b');
model = changeRxnBounds(model,'R_M_biomass_LLA_atpm',0,'b');
model = changeRxnBounds(model,'R_M_PROTS_LLA',0,'b');
model = changeRxnBounds(model,'R_M_PROTS_LLA_v2',0,'b');
model = changeRxnBounds(model,'R_M_PROTS_LLA_v3',0,'b');
model = changeRxnBounds(model,'R_M_MGt2pp_rvs',0,'b');%block infinite h[e]

% Block other glucose transporters
model = changeRxnBounds(model,'R_M_GLCpts_2',0,'b');
model = changeRxnBounds(model,'R_M_GLCpermease_fwd',0,'b');

%% Loop for dilution rate of 0.15 0.3 0.45 0.5 and 0.6.
[~, ~, exchange_raw] = xlsread('Exchange_reaction_setting.xlsx','Exp_bounds');

Exchange_reactions = exchange_raw(2:end,1);
header = exchange_raw(1,1:21);

mu_list = [0.15 0.3 0.45 0.5 0.6];

% without saturation factor
fluxes_global_saturation_factor_unchanged = [];
for i = 1:length(mu_list)
    mu = mu_list(i);
    factor_k = 1;
    % Set bounds for some exchange reactions.
    idx = find(contains(header,num2str(mu)));
    replicate = length(idx)/2;
    
	for j = 1:replicate
        lb_idx = idx(2*j-1);
        ub_idx = idx(2*j);
        LB = cell2mat(exchange_raw(2:end,lb_idx));
        UB = cell2mat(exchange_raw(2:end,ub_idx));
        model_tmp = changeRxnBounds(model,'R_biomass_dilution',mu,'b');
        model_tmp = changeRxnBounds(model_tmp,Exchange_reactions,LB,'l');
        model_tmp = changeRxnBounds(model_tmp,Exchange_reactions,UB,'u');
        disp(['Without sf: mu = ' num2str(mu) '; replicate = ' num2str(j)]);
        fileName = WriteLPSatFactor(model_tmp,mu,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,factor_k,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);

        command =sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,~,solME_full] = ReadSoplexResult(fileName_out,model_tmp);
        fluxes_global_saturation_factor_unchanged = [fluxes_global_saturation_factor_unchanged solME_full];
	end
end

% with saturation factor
load('Egsf1_result.mat');
mulist = global_saturation_factor_list(:,1);
sflist = global_saturation_factor_list(:,2);

fluxes_global_saturation_factor_changed = [];
for i = 1:length(mu_list)
    mu = mu_list(i);
    % Set bounds for some exchange reactions.
    idx = find(contains(header,num2str(mu)));
    replicate = length(idx)/2;
    sflisttmp = sflist(mulist == mu);
	for j = 1:replicate
        lb_idx = idx(2*j-1);
        ub_idx = idx(2*j);
        LB = cell2mat(exchange_raw(2:end,lb_idx));
        UB = cell2mat(exchange_raw(2:end,ub_idx));
        model_tmp = changeRxnBounds(model,'R_biomass_dilution',mu,'b');
        model_tmp = changeRxnBounds(model_tmp,Exchange_reactions,LB,'l');
        model_tmp = changeRxnBounds(model_tmp,Exchange_reactions,UB,'u');
        
        factor_k = sflisttmp(j);
        disp(['With sf: mu = ' num2str(mu) '; replicate = ' num2str(j)]);
        fileName = WriteLPSatFactor(model_tmp,mu,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,factor_k,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);

        command =sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,~,solME_full] = ReadSoplexResult(fileName_out,model_tmp);
        fluxes_global_saturation_factor_changed = [fluxes_global_saturation_factor_changed solME_full];
	end
end

cd Results/;
save('Cfd1_fluxes_without_sf.mat','fluxes_global_saturation_factor_unchanged');
save('Cfd1_fluxes_with_sf.mat','fluxes_global_saturation_factor_changed');
cd ../;

toc;