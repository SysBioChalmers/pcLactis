%% Estimate fraction of glucose transporter in total proteome

% Timing: ~ 4000 s

% By changing the fraction of glucose transporter in total proteome, we can
% study how glucose transporter affects the maximal growth rate.

% Simulated results will be saved in the folder 'Results'.

tic;
load('pcLactis_Model.mat');
model = pcLactis_Model;

%% Optimization setting.

rxnID = 'R_dummy_assumed_Monomer';
osenseStr = 'Maximize';

%% Parameters.
GAM = 40; %ATP coefficient in the new biomass equation.
NGAM = 3.5; %(mmol/gCDW/h)
f_unmodeled = 0.4; %proportion of unmodeled protein in total protein (g/g)

model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,f] = ChangeUnmodeledProtein(model,f_unmodeled);

kcat_glc = 180;%kcat value of glucose transporter (PMID: 28755958)
factor_k = 1;%global saturation factor

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

%% Main simulations.

% f_transporter_range = [0.002:0.002:0.01,0.1,0.2];
f_transporter_range = 0.5;
res = zeros(length(f_transporter_range),8);

for i = 1:length(f_transporter_range)
    f_transporter = f_transporter_range(i);
    
    mu_low = 0;
    mu_high = 0.8;
    
    while mu_high-mu_low > 0.00001
        mu_mid = (mu_low+mu_high)/2;
        disp(['f_transporter = ' num2str(f_transporter) '; mu = ' num2str(mu_mid)]);
        model = changeRxnBounds(model,'R_biomass_dilution',mu_mid,'b');
        model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*mu_mid,'l');
        
        fileName = WriteLP(model,mu_mid,f,osenseStr,rxnID,factor_k,...
                           f_transporter,kcat_glc,...
                           Info_enzyme,...
                           Info_mRNA,...
                           Info_protein,...
                           Info_ribosome,...
                           Info_tRNA);

        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t500 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
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
	
	fileName = WriteLP(model,mu_low,f,osenseStr,rxnID,factor_k,...
                       f_transporter,kcat_glc,...
                       Info_enzyme,...
                       Info_mRNA,...
                       Info_protein,...
                       Info_ribosome,...
                       Info_tRNA);
                   
	command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t500 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
    system(command,'-echo');
    fileName_out = 'Simulation.lp.out';
    [~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model);
    
    if strcmp(solME_status,'optimal')
        glc = -solME_full(strcmp(model.rxns,'R_M_EX_glc__D_e'));
        mu = solME_full(strcmp(model.rxns,'R_biomass_dilution'));
        ac = solME_full(strcmp(model.rxns,'R_M_EX_ac_e'));
        eth = solME_full(strcmp(model.rxns,'R_M_EX_etoh_e'));
        form = solME_full(strcmp(model.rxns,'R_M_EX_for_e'));
        lac = solME_full(strcmp(model.rxns,'R_M_EX_lac__L_e'));
        glc_transporter_weight = sum(CalculateEnzymeWeight(model,...
                                    {'M_GLCpts_1_Enzyme_c';...
                                     'M_GLCpts_2_Enzyme_c';...
                                     'M_GLCt2_fwd_Enzyme_c';...
                                     'M_GLCt2_rvs_Enzyme_c'},...
                                     solME_full,Info_enzyme));
        
        res(i,:) = [glc ac eth form lac mu glc_transporter_weight f_transporter];
    else
        res(i,:) = [0 0 0 0 0 0 0 f_transporter];
    end
end

estimated_f_transporter = max(res(:,7)) / 0.46;

clear ac ans command eth Exchange_AAs f f_transporter f_transporter_range;
clear f_unmodeled factor_k fileName fileName_out form GAM glc osenseStr;
clear glc_transporter_weight i kcat_glc lac LBfactor_AAs pcLactis_Model mu rxnID;
clear mu_high mu_low mu_mid NGAM solME_full solME_status;
clear Info_enzyme Info_mRNA Info_protein Info_protein Info_ribosome Info_tRNA;

cd Results/;
save('Egt_result.mat','res');
cd ../;

%% Figure
load('Egt_result.mat');
figure('Name','1');
hold on;
box on;
% x = res(:,8);
% y = res(:,6);
x = res(1:end-1,8);
y = res(1:end-1,6);
plot(x,y,'-o','LineWidth',1,'MarkerSize',1);
set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Fraction of glucose transporter','FontSize',14,'FontName','Helvetica');
ylabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');

% xlim([0 0.2]);
ylim([0 0.8]);

set(gcf,'position',[200 0 230 130]);
set(gca,'position',[0.23 0.25 0.72 0.7]);

clear x y;
toc;