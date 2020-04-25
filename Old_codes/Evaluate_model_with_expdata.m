%% Evaluate model using experimental data

% Timing: ~ 400 s

% In order to evaluate the model with the determined parameters, we used
% physiological data under chemostat conditions (PMID: 25828364) as
% constraints for exchange reactions to see if feasible solutions can be
% obtained.

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
factor_k = 1;%global saturation factor
f_transporter = 0.0082;%fraction of glucose transporter in total proteome

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

tf_res = [];
fluxes = zeros(length(model.rxns),0);
for i = 1:length(mu_list)
    mu = mu_list(i);
    
    model_tmp = changeRxnBounds(model,'R_biomass_dilution',mu,'b');
    
    % Set bounds for glucose transporter.
    glcT = f_transporter * 0.46 * mu / 1.0500867e2; % 1.0500867e2 is MW of glucose transporter in g/mmol
    model_tmp = changeRxnBounds(model_tmp,'R_dilution_M_GLCpts_2_Enzyme',glcT*0.999999,'b');
    
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
        
        fileName = WriteLPSatFactor(model_tmp,mu,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,1,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model_tmp);
        
        if strcmp(solME_status,'optimal')
            fluxes = cat(2,fluxes,solME_full);
        else
            fluxes = cat(2,fluxes,zeros(length(model_tmp.rxns),1));
        end
        
        tf_res = [tf_res;strcmp(solME_status,'optimal')];
        
	end
end
fluxes_global_saturation_factor_unchanged = fluxes;
cd Results/;
save('Emwe_fluxes.mat','fluxes_global_saturation_factor_unchanged');
cd ../;

flux_res = fluxes_global_saturation_factor_unchanged;
[~, excel_input, ~] = xlsread('Allocation.xlsx');
[~, n] = size(flux_res);

unmodeled_weight = 0.46 * f_unmodeled;

allocation_value = zeros(5,n);
total_proteome = zeros(1,n);
total_RNA = zeros(1,n);
total_rProtein = zeros(1,n);
total_inactive_enzyme = zeros(1,n);
mu_list = zeros(1,n);
total_Ribo = zeros(1,n);
q_glc = zeros(1,n);
total_rRNA = zeros(1,n);
total_tRNA = zeros(1,n);
total_mRNA = zeros(1,n);

for i = 1:n
    
    disp([num2str(i) '/' num2str(n)]);
    
    solME_full = flux_res(:,i);
    
    mu = solME_full(strcmp(model.rxns,'R_biomass_dilution'));
    mu_list(1,i) = mu;
    q_glc(1,i) = abs(solME_full(strcmp(model.rxns,'R_M_EX_glc__D_e'),:));
    
    [pathway, absvalue] = CalculateProteinAllocation(model,solME_full,Info_enzyme,excel_input);
    allocation_value(:,i) = absvalue;
    
    [Protein,RNA,rProtein,rRNA,tRNA,mRNA] = CalculateProteinAndRNA(model,solME_full,...
                                                             Info_enzyme,...
                                                             Info_ribosome,...
                                                             Info_tRNA,...
                                                             Info_mRNA);
    total_proteome(1,i) = Protein;
    total_RNA(1,i) = RNA;
    total_rProtein(1,i) = rProtein;
    total_Ribo(1,i) = rProtein + rRNA;
    
    [~,~,f_enzyme_inact,~] = CheckInactiveEnzyme(model,solME_full,factor_k);
    total_inactive_enzyme(1,i) = f_enzyme_inact;
    total_rRNA(1,i) = rRNA;
    total_tRNA(1,i) = tRNA;
    total_mRNA(1,i) = mRNA;
end
allocation_pathway = [pathway;'Ribosomal protein';'Inactive enzyme'];
allocation_value = [allocation_value;total_rProtein;total_inactive_enzyme];

allocation_pathway = [allocation_pathway;'Other modeled protein'];
allocation_value = [allocation_value;total_proteome-sum(allocation_value)-unmodeled_weight];

allocation_pathway = [allocation_pathway;'Unmodeled protein'];
allocation_value = [allocation_value;unmodeled_weight*ones(1,n)];

allocation_percentage = allocation_value./total_proteome;

toc;