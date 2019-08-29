load('Sglc2_fluxes_with_sf.mat');


% Experimental data
[~, ~, exp_raw] = xlsread('Exchange_reaction_setting.xlsx','Exchange_rates');
exp_mu = cell2mat(exp_raw(2,14:18));

exp_glc = -1*cell2mat(exp_raw(3,14:18));
exp_ac = cell2mat(exp_raw(4,14:18));
exp_eth = cell2mat(exp_raw(5,14:18));
exp_form = cell2mat(exp_raw(6,14:18));
exp_lac = cell2mat(exp_raw(7,14:18));

% Simulated data
load('pcLactis_Model.mat');
model = pcLactis_Model;

mu2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_biomass_dilution'),:);
glc2 = -fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_glc_LPAREN_e_RPAREN_'),:);
ac2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_ac_LPAREN_e_RPAREN_'),:);
eth2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_etoh_LPAREN_e_RPAREN_'),:);
form2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_for_LPAREN_e_RPAREN_'),:);
lac2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_lac_L_LPAREN_e_RPAREN_'),:);


% Figures
color_glc = [228,26,28]/255;
color_ac = [55,126,184]/255;
color_eth = [255,127,0]/255;
color_form = [77,175,74]/255;
color_lac = [152,78,163]/255;

% mixed to lactic acid with sf
figure('Name','3');
subplot(5,1,1);
hold on;
box on;
plot(mu2,glc2,'-','LineWidth',0.75,'Color',color_glc);
plot(exp_mu,exp_glc,'-.','LineWidth',0.75,'Color',color_glc);
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Flux','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('Glucose','FontSize',14,'FontName','Helvetica');

subplot(5,1,2);
hold on;
box on;
plot(mu2,ac2,'-','LineWidth',0.75,'Color',color_ac);
plot(exp_mu,exp_ac,'-.','LineWidth',0.75,'Color',color_ac);
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Flux','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('Acetate','FontSize',14,'FontName','Helvetica');

subplot(5,1,3);
hold on;
box on;
plot(mu2,eth2,'-','LineWidth',0.75,'Color',color_eth);
plot(exp_mu,exp_eth,'-.','LineWidth',0.75,'Color',color_eth);
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Flux','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('Ethanol','FontSize',14,'FontName','Helvetica');

subplot(5,1,4);
hold on;
box on;
plot(mu2,form2,'-','LineWidth',0.75,'Color',color_form);
plot(exp_mu,exp_form,'-.','LineWidth',0.75,'Color',color_form);
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Flux','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('Formate','FontSize',14,'FontName','Helvetica');

subplot(5,1,5);
hold on;
box on;
plot(mu2,lac2,'-','LineWidth',0.75,'Color',color_lac);
plot(exp_mu,exp_lac,'-.','LineWidth',0.75,'Color',color_lac);
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Flux','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('Lactate','FontSize',14,'FontName','Helvetica');