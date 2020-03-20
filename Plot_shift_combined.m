

load('Sglc_fluxes.mat');
flux_res_1 = fluxes_simulated_without_sf;

load('pcLactis_Model.mat');
model = pcLactis_Model;
GAM = 36; %ATP coefficient in the new biomass equation.
NGAM = 3; %(mmol/gCDW/h)
f_unmodeled = 0.4; %proportion of unmodeled protein in total protein (g/g)
model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,~] = ChangeUnmodeledProtein(model,f_unmodeled);
unmodeled_weight = 0.46 * f_unmodeled;
clear pcLactis_Model GAM NGAM f_unmodeled;

%% Metabolic shift
% Experimental data
[~, ~, exp_raw] = xlsread('Exchange_reaction_setting.xlsx','Exchange_rates');
exp_mu = cell2mat(exp_raw(2,14:18));

exp_glc = cell2mat(exp_raw(3,14:18));
exp_ac = cell2mat(exp_raw(4,14:18));
exp_eth = cell2mat(exp_raw(5,14:18));
exp_form = cell2mat(exp_raw(6,14:18));
exp_lac = cell2mat(exp_raw(7,14:18));
sd_glc = cell2mat(exp_raw(3,19:23));
sd_ac = cell2mat(exp_raw(4,19:23));
sd_eth = cell2mat(exp_raw(5,19:23));
sd_form = cell2mat(exp_raw(6,19:23));
sd_lac = cell2mat(exp_raw(7,19:23));

exp_arg = cell2mat(exp_raw(9,14:18));
exp_orn = cell2mat(exp_raw(8,14:18));
exp_nh4 = cell2mat(exp_raw(10,14:18));
sd_orn = cell2mat(exp_raw(8,19:23));
sd_arg = cell2mat(exp_raw(9,19:23));
sd_nh4 = cell2mat(exp_raw(10,19:23));

mu1 = flux_res_1(strcmp(model.rxns,'R_biomass_dilution'),:);
glc1 = flux_res_1(strcmp(model.rxns,'R_M_EX_glc__D_e'),:);
ac1 = flux_res_1(strcmp(model.rxns,'R_M_EX_ac_e'),:);
eth1 = flux_res_1(strcmp(model.rxns,'R_M_EX_etoh_e'),:);
form1 = flux_res_1(strcmp(model.rxns,'R_M_EX_for_e'),:);
lac1 = flux_res_1(strcmp(model.rxns,'R_M_EX_lac__L_e'),:);
orn1 = flux_res_1(strcmp(model.rxns,'R_M_EX_orn_e'),:);
nh41 = flux_res_1(strcmp(model.rxns,'R_M_EX_nh4_e'),:);
arg1 = flux_res_1(strcmp(model.rxns,'R_M_EX_arg__L_e'),:);

% Figures

color_glc = [247,129,191]/255;
color_ac = [55,126,184]/255;
% color_eth = [255,127,0]/255;
% color_form = [77,175,74]/255;
color_eth = color_ac;
color_form = color_ac;
color_lac = [228,26,28]/255;
% color_arg = [152,78,163]/255;
color_arg = [0,0,0];

figure('Name','1');
subplot(6,1,1);
hold on;
box on;
plot(mu1,glc1,'-','LineWidth',0.5,'Color',color_glc);
scatter(mu1,glc1,5,'o','LineWidth',1,'MarkerEdgeColor',color_glc,'MarkerEdgeAlpha',0.8);
scatter(exp_mu,exp_glc,10,'x','LineWidth',1,'MarkerEdgeColor',color_glc);
x = exp_mu; y = exp_glc; yu = exp_glc + sd_glc; yl = exp_glc - sd_glc;
plot(x,y,':','LineWidth',1,'Color',color_glc);
fill([x fliplr(x)],[yu fliplr(yl)],color_glc,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
ylim([-25 0]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Flux','FontSize',7,'FontName','Helvetica');
title('Glucose','FontSize',6,'FontName','Helvetica');

subplot(6,1,2);
hold on;
box on;
plot(mu1,ac1,'-','LineWidth',0.5,'Color',color_ac);
scatter(mu1,ac1,5,'o','LineWidth',1,'MarkerEdgeColor',color_ac,'MarkerEdgeAlpha',0.8);
scatter(exp_mu,exp_ac,10,'x','LineWidth',1,'MarkerEdgeColor',color_ac);
x = exp_mu; y = exp_ac; yu = exp_ac + sd_ac; yl = exp_ac - sd_ac;
plot(x,y,':','LineWidth',1,'Color',color_ac);
fill([x fliplr(x)],[yu fliplr(yl)],color_ac,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
ylim([0 13]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Flux','FontSize',7,'FontName','Helvetica');
title('Acetate','FontSize',6,'FontName','Helvetica');

subplot(6,1,3);
hold on;
box on;
plot(mu1,eth1,'-','LineWidth',0.5,'Color',color_eth);
scatter(mu1,eth1,5,'o','LineWidth',1,'MarkerEdgeColor',color_eth,'MarkerEdgeAlpha',0.8);
scatter(exp_mu,exp_eth,10,'x','LineWidth',1,'MarkerEdgeColor',color_eth);
x = exp_mu; y = exp_eth; yu = exp_eth + sd_eth; yl = exp_eth - sd_eth;
plot(x,y,':','LineWidth',1,'Color',color_eth);
fill([x fliplr(x)],[yu fliplr(yl)],color_eth,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
ylim([0 13]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Flux','FontSize',7,'FontName','Helvetica');
title('Ethanol','FontSize',6,'FontName','Helvetica');

subplot(6,1,4);
hold on;
box on;
plot(mu1,form1,'-','LineWidth',0.5,'Color',color_form);
scatter(mu1,form1,5,'o','LineWidth',1,'MarkerEdgeColor',color_form,'MarkerEdgeAlpha',0.8);
scatter(exp_mu,exp_form,10,'x','LineWidth',1,'MarkerEdgeColor',color_form);
x = exp_mu; y = exp_form; yu = exp_form + sd_form; yl = exp_form - sd_form;
plot(x,y,':','LineWidth',1,'Color',color_form);
fill([x fliplr(x)],[yu fliplr(yl)],color_form,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
ylim([0 25]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Flux','FontSize',7,'FontName','Helvetica');
title('Formate','FontSize',6,'FontName','Helvetica');

subplot(6,1,5);
hold on;
box on;
plot(mu1,lac1,'-','LineWidth',0.5,'Color',color_lac);
scatter(mu1,lac1,5,'o','LineWidth',1,'MarkerEdgeColor',color_lac,'MarkerEdgeAlpha',0.8);
scatter(exp_mu,exp_lac,10,'x','LineWidth',1,'MarkerEdgeColor',color_lac);
x = exp_mu; y = exp_lac; yu = exp_lac + sd_lac; yl = exp_lac - sd_lac;
plot(x,y,':','LineWidth',1,'Color',color_lac);
fill([x fliplr(x)],[yu fliplr(yl)],color_lac,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
ylim([0 50]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Flux','FontSize',7,'FontName','Helvetica');
title('Lactate','FontSize',6,'FontName','Helvetica');

subplot(6,1,6);
hold on;
box on;
plot(mu1,arg1,'-','LineWidth',0.5,'Color',color_arg);
scatter(mu1,arg1,5,'o','LineWidth',1,'MarkerEdgeColor',color_arg,'MarkerEdgeAlpha',0.8);
scatter(exp_mu,exp_arg,10,'x','LineWidth',1,'MarkerEdgeColor',color_arg);
x = exp_mu; y = exp_arg; yu = exp_arg + sd_arg; yl = exp_arg - sd_arg;
plot(x,y,':','LineWidth',1,'Color',color_arg);
fill([x fliplr(x)],[yu fliplr(yl)],color_arg,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
ylim([-1.3 0]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Flux','FontSize',7,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
title('Arginine','FontSize',6,'FontName','Helvetica');
legend();

set(gcf,'position',[500 270 90 300]);
clear;