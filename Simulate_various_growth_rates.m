%% Simulate glucose-limited chemostats (objective: minimizing glucose uptake rate)

% Timing: ~  s

% Without and with the saturation saturation factor are performed.

% Simulated results will be saved in the folder 'Results'.

% Figures can be obtained by running the codes starting from line 157.

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
f_transporter = 1;%no constraint on glucose transporter

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

%% Main part.

D_list = 0.1:0.1:0.7;%unit: /h

% without saturation factor
factor_k = 1;%global saturation factor
factor_glc = 1;%assume glucose transporter is fully saturated

fluxes_simulated_without_sf = zeros(length(model.rxns),length(D_list));

for i = 1:length(D_list)
    
    D = D_list(i);
	disp(['Without sf: D = ' num2str(D)]);
	model = changeRxnBounds(model,'R_biomass_dilution',D,'b');
	model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*D,'l');
	fileName = WriteLPSatFactor(model,D,f,osenseStr,rxnID,factor_k,...
                                f_transporter,kcat_glc,factor_glc,...
                                Info_enzyme,...
                                Info_mRNA,...
                                Info_protein,...
                                Info_ribosome,...
                                Info_tRNA);
                   
	command =sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
	system(command,'-echo');
	fileName_out = 'Simulation.lp.out';
	[~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model);
    
    if strcmp(solME_status,'optimal')
        fluxes_simulated_without_sf(:,i) = solME_full;
    else
        fluxes_simulated_without_sf(:,i) = zeros(length(model.rxns),1);
    end
end

% with saturation factor
% obtain the global saturation factor
load('Egsf1_result.mat');
x = global_saturation_factor_list(:,1);
y = global_saturation_factor_list(:,2);
x = x(y ~= 1);
y = y(y ~= 1);
sf_coeff = x\y;
clear x y;

fluxes_simulated_with_sf = zeros(length(model.rxns),length(D_list));

for i = 1:length(D_list)
    
    D = D_list(i);
    
	factor_k = sf_coeff * D;
    if factor_k > 1
        factor_k = 1;
    end
    factor_glc = factor_k;
    
	disp(['With sf: D = ' num2str(D)]);
	model = changeRxnBounds(model,'R_biomass_dilution',D,'b');
	model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*D,'l');
    
	fileName = WriteLPSatFactor(model,D,f,osenseStr,rxnID,factor_k,...
                                f_transporter,kcat_glc,factor_glc,...
                                Info_enzyme,...
                                Info_mRNA,...
                                Info_protein,...
                                Info_ribosome,...
                                Info_tRNA);
                   
	command =sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
	system(command,'-echo');
	fileName_out = 'Simulation.lp.out';
	[~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model);
    
    if strcmp(solME_status,'optimal')
        fluxes_simulated_with_sf(:,i) = solME_full;
    else
        fluxes_simulated_with_sf(:,i) = zeros(length(model.rxns,1));
    end
end

cd Results/;
save('Svgr_fluxes_without_sf.mat','fluxes_simulated_without_sf');
save('Svgr_fluxes_with_sf.mat','fluxes_simulated_with_sf');
cd ../;

clear;
toc;

%% Figure.
load('Svgr_fluxes_without_sf.mat');
load('Svgr_fluxes_with_sf.mat');

%Metabolic shift
% Experimental data
[~, ~, exp_raw] = xlsread('Exchange_reaction_setting.xlsx','Exchange_rates');
exp_mu = cell2mat(exp_raw(2,12:16));

f_exp_mix = cell2mat(exp_raw(15,12:16));
f_exp_lac = cell2mat(exp_raw(14,12:16));
f_exp_mix_sd = cell2mat(exp_raw(15,17:21));
f_exp_lac_sd = cell2mat(exp_raw(14,17:21));

% Simulated data
load('pcLactis_Model.mat');
model = pcLactis_Model;
mu1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_biomass_dilution'),:);
glc1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_glc_LPAREN_e_RPAREN_'),:);
ac1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_ac_LPAREN_e_RPAREN_'),:);
eth1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_etoh_LPAREN_e_RPAREN_'),:);
form1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_for_LPAREN_e_RPAREN_'),:);
lac1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_lac_L_LPAREN_e_RPAREN_'),:);
f_mix1 = round((ac1*2+eth1*2+form1)./(-glc1*6),2);
f_lac1 = (lac1*3)./(-glc1*6);

mu2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_biomass_dilution'),:);
glc2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_glc_LPAREN_e_RPAREN_'),:);
ac2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_ac_LPAREN_e_RPAREN_'),:);
eth2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_etoh_LPAREN_e_RPAREN_'),:);
form2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_for_LPAREN_e_RPAREN_'),:);
lac2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_lac_L_LPAREN_e_RPAREN_'),:);
f_mix2 = (ac2*2+eth2*2+form2)./(-glc2*6);
f_lac2 = (lac2*3)./(-glc2*6);

% Figures
% mixed to lactic acid without sf
color_mixed = [94,60,153]/255;
color_lactic = [230,97,1]/255;
figure('Name','1');
hold on;
box on;
errorbar(exp_mu,f_exp_lac,f_exp_lac_sd,'--o','LineWidth',0.5,'Color',color_lactic,'MarkerSize',7,'MarkerEdgeColor',color_lactic);
errorbar(exp_mu,f_exp_mix,f_exp_mix_sd,'--o','LineWidth',0.5,'Color',color_mixed,'MarkerSize',7,'MarkerEdgeColor',color_mixed);
plot(mu1,f_lac1,'-o','LineWidth',0.5,'Color',color_lactic,'MarkerSize',7,'MarkerEdgeColor',color_lactic,'MarkerFaceColor',color_lactic);
plot(mu1,f_mix1,'-o','LineWidth',0.5,'Color',color_mixed,'MarkerSize',7,'MarkerEdgeColor',color_mixed,'MarkerFaceColor',color_mixed);
set(gca,'YTick',0:0.2:1);
set(gca,'ycolor','k');
set(gca,'XTick',0.1:0.2:0.7);
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Fraction in glucose','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
legend({'Sim mixed acid',...
        'Sim lactate'...
        'Exp mixed acid',...
        'Exp lactate'},'FontSize',12,'FontName','Helvetica','location','north');
set(gcf,'position',[500 200 200 180]);
set(gca,'position',[0.21 0.23 0.73 0.75]);

% mixed to lactic acid with sf
figure('Name','2');
hold on;
box on;
errorbar(exp_mu,f_exp_lac,f_exp_lac_sd,'--o','LineWidth',0.5,'Color',color_lactic,'MarkerSize',7,'MarkerEdgeColor',color_lactic);
errorbar(exp_mu,f_exp_mix,f_exp_mix_sd,'--o','LineWidth',0.5,'Color',color_mixed,'MarkerSize',7,'MarkerEdgeColor',color_mixed);
plot(mu2,f_lac2,'-o','LineWidth',0.5,'Color',color_lactic,'MarkerSize',7,'MarkerEdgeColor',color_lactic,'MarkerFaceColor',color_lactic);
plot(mu2,f_mix2,'-o','LineWidth',0.5,'Color',color_mixed,'MarkerSize',7,'MarkerEdgeColor',color_mixed,'MarkerFaceColor',color_mixed);
set(gca,'YTick',0:0.2:1);
set(gca,'ycolor','k');
set(gca,'XTick',0.1:0.2:0.7);
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Fraction in glucose','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
% legend({'Sim mixed acid',...
%         'Sim lactate'...
%         'Exp mixed acid',...
%         'Exp lactate'},'FontSize',12,'FontName','Helvetica','location','north');
set(gcf,'position',[500 0 200 180]);
set(gca,'position',[0.21 0.23 0.73 0.75]);

clear;
