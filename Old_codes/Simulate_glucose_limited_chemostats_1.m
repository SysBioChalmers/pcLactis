%% Simulate glucose-limited chemostats (objective: minimizing glucose uptake rate)

% Timing: ~ 1200 s

% Without and with the saturation saturation factor are performed.

% Simulated results will be saved in the folder 'Results'.

% Figures can be obtained by running the codes starting from line 151.

tic;
load('pcLactis_Model.mat');
model = pcLactis_Model;

%% Optimization setting.

rxnID = 'R_M_EX_glc__D_e';
osenseStr = 'Maximize';

%% Parameters.
GAM = 36;%ATP coefficient in the new biomass equation.
NGAM = 3; %(mmol/gCDW/h)
f_unmodeled = 0.4; %proportion of unmodeled protein in total protein (g/g)

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
model = changeRxnBounds(model,'R_M_PROTS_LLA_v3',0,'b');

% Block other glucose transporters
model = changeRxnBounds(model,'R_M_GLCpts_1',0,'b');
model = changeRxnBounds(model,'R_M_GLCt2_fwd',0,'b');

% Block one of ADH isozymes llmg_0955
model = changeRxnBounds(model,'R_M_ALCD2x_1_rvs',0,'b');

%% Main part.

D_list = 0.1:0.1:0.7;%unit: /h

% without saturation factor
factor_k = 1;%global saturation factor
fluxes_simulated_without_sf = zeros(length(model.rxns),length(D_list));

for i = 1:length(D_list)
    
    D = D_list(i);
    disp(['Without sf: D = ' num2str(D)]);
	model = changeRxnBounds(model,'R_biomass_dilution',D,'b');
	model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*D,'l');
	fileName = WriteLPSatFactor(model,D,f,osenseStr,rxnID,factor_k,...
                                f_transporter,kcat_glc,factor_k,...
                                Info_enzyme,...
                                Info_mRNA,...
                                Info_protein,...
                                Info_ribosome,...
                                Info_tRNA);
                   
	command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
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

	disp(['With sf: D = ' num2str(D)]);
	model = changeRxnBounds(model,'R_biomass_dilution',D,'b');
	model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*D,'l');
	fileName = WriteLPSatFactor(model,D,f,osenseStr,rxnID,factor_k,...
                                f_transporter,kcat_glc,factor_k,...
                                Info_enzyme,...
                                Info_mRNA,...
                                Info_protein,...
                                Info_ribosome,...
                                Info_tRNA);
                   
	command =sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
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
save('Sglc1_fluxes_without_sf.mat','fluxes_simulated_without_sf');
save('Sglc1_fluxes_with_sf.mat','fluxes_simulated_with_sf');
cd ../;

clear;
toc;

%% Figure.
load('Sglc1_fluxes_without_sf.mat');
load('Sglc1_fluxes_with_sf.mat');

%Metabolic shift
% Experimental data
[~, ~, exp_raw] = xlsread('Exchange_reaction_setting.xlsx','Exchange_rates');
exp_mu = cell2mat(exp_raw(2,14:18));

f_exp_mix = cell2mat(exp_raw(15,14:18));
f_exp_lac = cell2mat(exp_raw(14,14:18));
f_exp_mix_sd = cell2mat(exp_raw(15,19:23));
f_exp_lac_sd = cell2mat(exp_raw(14,19:23));

% Simulated data
load('pcLactis_Model.mat');
model = pcLactis_Model;
mu1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_biomass_dilution'),:);
glc1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_glc__D_e'),:);
ac1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_ac_e'),:);
eth1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_etoh_e'),:);
form1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_for_e'),:);
lac1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_lac__L_e'),:);
f_mix1 = (ac1*2+eth1*2+form1)./(-glc1*6);
f_lac1 = (lac1*3)./(-glc1*6);

mu2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_biomass_dilution'),:);
glc2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_glc__D_e'),:);
ac2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_ac_e'),:);
eth2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_etoh_e'),:);
form2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_for_e'),:);
lac2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_lac__L_e'),:);
f_mix2 = (ac2*2+eth2*2+form2)./(-glc2*6);
f_lac2 = (lac2*3)./(-glc2*6);

% Figures
% mixed to lactic acid without sf
color_mixed = [94,60,153]/255;
color_lactic = [230,97,1]/255;
figure('Name','2');
hold on;
box on;

plot(mu1,f_lac1,'-','LineWidth',0.75,'Color',color_lactic);
x = exp_mu; y = f_exp_lac; yu = f_exp_lac + f_exp_lac_sd; yl = f_exp_lac - f_exp_lac_sd;
plot(x,y,'-.','LineWidth',0.75,'Color',color_lactic);
fill([x fliplr(x)],[yu fliplr(yl)],color_lactic,'linestyle','none','FaceAlpha',0.3);

plot(mu1,f_mix1,'-','LineWidth',0.75,'Color',color_mixed);
x = exp_mu; y = f_exp_mix; yu = f_exp_mix + f_exp_mix_sd; yl = f_exp_mix - f_exp_mix_sd;
plot(x,y,'-.','LineWidth',0.75,'Color',color_mixed);
fill([x fliplr(x)],[yu fliplr(yl)],color_mixed,'linestyle','none','FaceAlpha',0.3);

% errorbar(exp_mu,f_exp_lac,f_exp_lac_sd,'-.','LineWidth',0.75,'Color',color_lactic);
% errorbar(exp_mu,f_exp_mix,f_exp_mix_sd,'-.','LineWidth',0.75,'Color',color_mixed);

set(gca,'YTick',0:0.2:1);
set(gca,'ycolor','k');
set(gca,'XTick',0.1:0.2:0.7);
xlim([0.1 0.7]);
ylim([-0.05 1.05]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Fraction in glucose','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
legend({'Sim lactic acid',...
        'Exp lactic acid'...
        'SD lactic acid'...
        'Sim mixed acid',...
        'Exp mixed acid'...
        'SD mixed acid'},'FontSize',12,'FontName','Helvetica','location','north');
set(gcf,'position',[500 200 200 180]);
set(gca,'position',[0.21 0.23 0.73 0.75]);

% mixed to lactic acid with sf
figure('Name','3');
hold on;
box on;

plot(mu2,f_lac2,'-','LineWidth',0.75,'Color',color_lactic);
x = exp_mu; y = f_exp_lac; yu = f_exp_lac + f_exp_lac_sd; yl = f_exp_lac - f_exp_lac_sd;
plot(x,y,'-.','LineWidth',0.75,'Color',color_lactic);
fill([x fliplr(x)],[yu fliplr(yl)],color_lactic,'linestyle','none','FaceAlpha',0.3);

plot(mu2,f_mix2,'-','LineWidth',0.75,'Color',color_mixed);
x = exp_mu; y = f_exp_mix; yu = f_exp_mix + f_exp_mix_sd; yl = f_exp_mix - f_exp_mix_sd;
plot(x,y,'-.','LineWidth',0.75,'Color',color_mixed);
fill([x fliplr(x)],[yu fliplr(yl)],color_mixed,'linestyle','none','FaceAlpha',0.3);

% errorbar(exp_mu,f_exp_lac,f_exp_lac_sd,'-.','LineWidth',0.75,'Color',color_lactic);
% errorbar(exp_mu,f_exp_mix,f_exp_mix_sd,'-.','LineWidth',0.75,'Color',color_mixed);

set(gca,'YTick',0:0.2:1);
set(gca,'ycolor','k');
set(gca,'XTick',0.1:0.2:0.7);
xlim([0.1 0.7]);
ylim([-0.05 1.05]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Fraction in glucose','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
set(gcf,'position',[500 0 200 180]);
set(gca,'position',[0.21 0.23 0.73 0.75]);

clear;
