
load('Sglc_fluxes.mat');
flux_res = fluxes_simulated_without_sf;

load('pcLactis_Model.mat');
model = pcLactis_Model;
factor_k = 1;

load('Sglc_result.mat');
Km = 21;%Km of glucose transporter, unit: uM (PMID: 30630406)
S = glc_conc_without_sf(:,2);
glcT_saturation = S./(S+Km);

%% Metabolic shift
% Experimental data
[~, ~, exp_raw] = xlsread('Exchange_reaction_setting.xlsx','Exchange_rates');
exp_mu = cell2mat(exp_raw(2,14:18));

exp_glc = -cell2mat(exp_raw(3,14:18));
exp_ac = cell2mat(exp_raw(4,14:18));
exp_form = cell2mat(exp_raw(6,14:18));
exp_lac = cell2mat(exp_raw(7,14:18));
sd_glc = cell2mat(exp_raw(3,19:23));
sd_ac = cell2mat(exp_raw(4,19:23));
sd_form = cell2mat(exp_raw(6,19:23));
sd_lac = cell2mat(exp_raw(7,19:23));

% Simulated data
mu = flux_res(strcmp(model.rxns,'R_biomass_dilution'),:);
glc = -flux_res(strcmp(model.rxns,'R_M_EX_glc__D_e'),:);
ac = flux_res(strcmp(model.rxns,'R_M_EX_ac_e'),:);
form = flux_res(strcmp(model.rxns,'R_M_EX_for_e'),:);
lac = flux_res(strcmp(model.rxns,'R_M_EX_lac__L_e'),:);
%inactive enzyme
[~,b] = size(flux_res);
total_inactive_enzyme = zeros(1,b);
for i = 1:b
    [~,~,f_enzyme_inact,~] = CheckInactiveEnzyme(model,flux_res(:,i),factor_k);
    if f_enzyme_inact > 1e-6
        total_inactive_enzyme(i) = f_enzyme_inact;
    end
end

% Figures
color_glc = [152,78,163]/255;
color_ac = [55,126,184]/255;
color_form = [77,175,74]/255;
color_lac = [228,26,28]/255;
color_glcT = [197,27,138]/255;
color_inactE = [197,27,138]/255;

figure('Name','1');
subplot(4,1,1);
hold on;
box on;
%glcuose
plot(mu,glc,'-','Marker','.','MarkerSize',5,'LineWidth',0.75,'Color',color_glc);
scatter(exp_mu,exp_glc,10,'o','LineWidth',1,'MarkerEdgeColor',color_glc);
% x = exp_mu; y = exp_glc; yu = exp_glc + sd_glc; yl = exp_glc - sd_glc;
% fill([x fliplr(x)],[yu fliplr(yl)],color_glc,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
% ylim([-25 0]);
set(gca, 'XColor','k');
set(gca, 'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
% ylabel(['Flux',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
title('Glucose uptake','FontSize',6,'FontName','Helvetica','Color','k');

subplot(4,1,2);
hold on;
box on;
%acetate
plot(mu,ac,'-','Marker','.','MarkerSize',5,'LineWidth',0.75,'Color',color_ac);
scatter(exp_mu,exp_ac,10,'o','LineWidth',1,'MarkerEdgeColor',color_ac);
% x = exp_mu; y = exp_ac; yu = exp_ac + sd_ac; yl = exp_ac - sd_ac;
% fill([x fliplr(x)],[yu fliplr(yl)],color_ac,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
% ylim([-25 0]);
set(gca, 'XColor','k');
set(gca, 'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
% ylabel(['Flux',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
title('Acetate production','FontSize',6,'FontName','Helvetica','Color','k');

subplot(4,1,3);
hold on;
box on;
%formate
plot(mu,form,'-','Marker','.','MarkerSize',5,'LineWidth',0.75,'Color',color_form);
scatter(exp_mu,exp_form,10,'o','LineWidth',1,'MarkerEdgeColor',color_form);
% x = exp_mu; y = exp_form; yu = exp_form + sd_form; yl = exp_form - sd_form;
% fill([x fliplr(x)],[yu fliplr(yl)],color_form,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
% ylim([-25 0]);
set(gca, 'XColor','k');
set(gca, 'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel(['Flux',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica','Color','k');
title('Formate production','FontSize',6,'FontName','Helvetica','Color','k');

subplot(4,1,4);
hold on;
box on;
%lactate
plot(mu,lac,'-','Marker','.','MarkerSize',5,'LineWidth',0.75,'Color',color_lac);
scatter(exp_mu,exp_lac,10,'o','LineWidth',1,'MarkerEdgeColor',color_lac);
% x = exp_mu; y = exp_lac; yu = exp_lac + sd_lac; yl = exp_lac - sd_lac;
% fill([x fliplr(x)],[yu fliplr(yl)],color_lac,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
% ylim([-25 0]);
set(gca, 'XColor','k');
set(gca, 'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
% ylabel(['Flux',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica','Color','k');
title('Lactate production','FontSize',6,'FontName','Helvetica','Color','k');
legend();
set(gcf,'position',[500 270 90 200]);

figure('Name','2');
subplot(1,3,1);
hold on;
box on;
%glucose transporter saturation
plot(mu,glcT_saturation,'-','Marker','.','MarkerSize',5,'LineWidth',0.75,'Color',color_glcT);
% exp_mumax = 0.71;
% exp_K = 0.019*1000;
% x = 0:1:3000;
% y = exp_mumax*x./(exp_K+x);
% plot(y,x./(x+Km),'-','LineWidth',0.5,'Color','k');
xlim([0.1 0.7]);
set(gca, 'XColor','k');
set(gca, 'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Saturation','FontSize',7,'FontName','Helvetica','Color','k');
title('Glucose transporter','FontSize',6,'FontName','Helvetica','Color','k');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica','Color','k');

subplot(1,3,2);
hold on;
box on;
%inavtive enzyme
plot(mu,total_inactive_enzyme,'-','Marker','.','MarkerSize',5,'LineWidth',0.75,'Color',color_inactE);
fill([[0.5 0.6] [0.6 0.5]],[[1 1] [0 0]],'k','linestyle','none','FaceAlpha',0.2);
xlim([0.1 0.7]);
ylim([0 0.2]);
set(gca, 'XColor','k');
set(gca, 'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel(['Content',char(13,10)','(g/gCDW)'],'FontSize',7,'FontName','Helvetica','Color','k');
title('Inactive enzyme','FontSize',6,'FontName','Helvetica','Color','k');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica','Color','k');

subplot(1,3,3);
hold on;
box on;
%inavtive enzyme
plot(mu,total_inactive_enzyme,'-','Marker','.','MarkerSize',10,'LineWidth',0.75,'Color',color_inactE);
fill([[0.5 0.6] [0.6 0.5]],[[1 1] [0 0]],'k','linestyle','none','FaceAlpha',0.2);
xlim([0.48 0.62]);
ylim([0 0.012]);
set(gca, 'XColor','k');
set(gca, 'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel(['Content',char(13,10)','(g/gCDW)'],'FontSize',7,'FontName','Helvetica','Color','k');
title('Inactive enzyme','FontSize',6,'FontName','Helvetica','Color','k');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica','Color','k');

set(gcf,'position',[200 300 290 65]);
