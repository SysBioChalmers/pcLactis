load('Sglc2_fluxes_with_sf.mat');
flux_res = fluxes_simulated_with_sf;
% load('Sglc2_fluxes_without_sf.mat');
% flux_res = fluxes_simulated_without_sf;

load('pcLactis_Model.mat');
model = pcLactis_Model;
GAM = 36;%ATP coefficient in the new biomass equation.
NGAM = 2; %(mmol/gCDW/h)
f_unmodeled = 0.4; %proportion of unmodeled protein in total protein (g/g)
model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,~] = ChangeUnmodeledProtein(model,f_unmodeled);
unmodeled_weight = 0.46 * f_unmodeled;
clear pcLactis_Model GAM NGAM f_unmodeled;

%% Data
mu = flux_res(strcmp(model.rxns,'R_biomass_dilution'),:);

% LDH
dil_ldh1 = flux_res(strcmp(model.rxns,'R_dilution_M_LDH_L_1_rvs_Enzyme'),:);
dil_ldh2 = flux_res(strcmp(model.rxns,'R_dilution_M_LDH_L_2_rvs_Enzyme'),:);
dil_ldh3 = flux_res(strcmp(model.rxns,'R_dilution_M_LDH_L_3_rvs_Enzyme'),:);
abs_ldh1 = dil_ldh1 ./ mu;
abs_ldh2 = dil_ldh2 ./ mu;
abs_ldh3 = dil_ldh3 ./ mu;

% ACK
dil_ack1 = flux_res(strcmp(model.rxns,'R_dilution_M_ACKr_1_rvs_Enzyme'),:);
dil_ack2 = flux_res(strcmp(model.rxns,'R_dilution_M_ACKr_2_rvs_Enzyme'),:);
abs_ack1 = dil_ack1 ./ mu;
abs_ack2 = dil_ack2 ./ mu;

% ADH
dil_adh1 = flux_res(strcmp(model.rxns,'R_dilution_M_ALCD2x_1_rvs_Enzyme'),:);
dil_adh2 = flux_res(strcmp(model.rxns,'R_dilution_M_ALCD2x_2_rvs_Enzyme'),:);
abs_adh1 = dil_adh1 ./ mu;
abs_adh2 = dil_adh2 ./ mu;

% PFL
dil_pfl = flux_res(strcmp(model.rxns,'R_dilution_M_PFL_fwd_Enzyme'),:);
abs_pfl = dil_pfl ./ mu;

% ArgA (arginine deiminase)
dil_arga = flux_res(strcmp(model.rxns,'R_dilution_M_ARGDI_Enzyme'),:);
abs_arga = dil_arga ./ mu;

%% Figures

figure('Name','1');
subplot(1,5,1);
hold on;
box on;
plot(mu,abs_ldh1,'k-o',mu,abs_ldh2,'k-o',mu,abs_ldh3,'k-o');
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Protein (mmol/gCDW)','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('LDH','FontSize',14,'FontName','Helvetica');

subplot(1,5,2);
hold on;
box on;
plot(mu,abs_ack1,'k-o',mu,abs_ack2,'k-o');
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Protein (mmol/gCDW)','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('ACK','FontSize',14,'FontName','Helvetica');

subplot(1,5,3);
hold on;
box on;
plot(mu,abs_adh1,'k-o',mu,abs_adh2,'k-o');
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Protein (mmol/gCDW)','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('ADH','FontSize',14,'FontName','Helvetica');

subplot(1,5,4);
hold on;
box on;
plot(mu,abs_pfl,'k-o');
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Protein (mmol/gCDW)','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('PFL','FontSize',14,'FontName','Helvetica');

subplot(1,5,5);
hold on;
box on;
plot(mu,abs_arga,'k-o');
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Protein (mmol/gCDW)','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('ArgA','FontSize',14,'FontName','Helvetica');

figure('Name','2');
mu_new = [0.3 0.4 0.45 0.5 0.6];
idx = [5 7 8 9 11];
% mu_new = [0.3 0.5 0.6];
% idx = [5 9 11];

subplot(1,3,1);
hold on;
box on;
plot(mu_new,log2(abs_ack2(idx)/abs_ack2(2)),'k-o');
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('log_2 Ratio','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('ACK','FontSize',14,'FontName','Helvetica');

subplot(1,3,2);
hold on;
box on;
plot(mu_new,log2(abs_adh2(idx)/abs_adh2(2)),'k-o');
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
% ylabel('log_2 Ratio','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('ADH','FontSize',14,'FontName','Helvetica');

subplot(1,3,3);
hold on;
box on;
plot(mu_new,log2(abs_pfl(idx)/abs_pfl(2)),'k-o');
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
% ylabel('log_2 Ratio','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('PFL','FontSize',14,'FontName','Helvetica');

% subplot(1,5,5);
% hold on;
% box on;
% plot(mu_new,log2(abs_arga(idx)/abs_arga(2)),'k-o');
% xlim([0.1 0.7]);
% set(gca,'FontSize',12,'FontName','Helvetica');
% ylabel('log_2 Ratio','FontSize',14,'FontName','Helvetica');
% xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
% title('ArgA','FontSize',14,'FontName','Helvetica');

set(gcf,'position',[400 400 500 120]);

