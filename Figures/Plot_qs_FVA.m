%% Figure.
load('Eguv_fluxes.mat');
load('pcLactis_Model.mat');
model = pcLactis_Model;

mu = fluxes_simulated(strcmp(model.rxns,'R_biomass_dilution'),:);
glc = -fluxes_simulated(strcmp(model.rxns,'R_M_EX_glc__D_e'),:);
mu_list = mu(1,1:2:end);
low_glc_list = glc(1,1:2:end);
high_glc_list = glc(1,2:2:end);
aver_glc_list = (low_glc_list + high_glc_list)/2;
rel_aver_list = 100 * ones(length(aver_glc_list),1);
rel_err_list = 100 * (high_glc_list - low_glc_list)/2 ./ aver_glc_list;

% Relative values
figure('Name','1');
errorbar(mu_list,rel_aver_list,rel_err_list,'.','LineWidth',1);

set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',8,'FontName','Helvetica');
ylabel('Relative q_s (%)','FontSize',8,'FontName','Helvetica');

ylim([99.9,100.1]);

set(gcf,'position',[0 300 230 130]);
set(gca,'position',[0.27 0.29 0.68 0.66]);

% Absolute values
figure('Name','2');
hold on;
box on;
plot(mu_list,aver_glc_list,'o','LineWidth',1);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',8,'FontName','Helvetica');
ylabel(['Absolute q_s',char(13,10)','(mmol/gCDW/h)'],'FontSize',8,'FontName','Helvetica');

set(gcf,'position',[0 0 230 130]);
set(gca,'position',[0.27 0.29 0.68 0.66]);

clear;