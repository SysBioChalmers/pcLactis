load('Sglc2_result_without_sf.mat');
load('Sglc2_result_with_sf.mat');
% glc_conc_without_sf = glc_conc_without_sf([1:2:7,8:1:length(glc_conc_without_sf)],:);
% glc_conc_with_sf = glc_conc_with_sf([1:2:7,8:1:length(glc_conc_with_sf)],:);
clr = [178,24,43]/255;

%Glucose concentration and growth rate
figure('Name','1');
hold on;
box on;

exp_mumax = 0.71;
exp_K = 0.019*1000;
x = 0:1:250;
y = exp_mumax*x./(exp_K+x);
plot(x,y,'k-','LineWidth',1);

y1 = glc_conc_without_sf(:,1);
x1 = glc_conc_without_sf(:,2);
scatter(x1,y1,20,'o','LineWidth',1,'MarkerEdgeColor',clr,'MarkerEdgeAlpha',0.8);
y2 = glc_conc_with_sf(:,1);
x2 = glc_conc_with_sf(:,2);
scatter(x2,y2,20,'o','filled','LineWidth',1,'MarkerEdgeColor',clr,'MarkerFaceColor',clr,'MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8);

set(gca,'FontSize',6,'FontName','Helvetica');
% legend({'Before','After'},'FontSize',7,'FontName','Helvetica','location','north');
xlabel('Glucose concentration (uM)','FontSize',7,'FontName','Helvetica');
ylabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');

xlim([0 240]);
ylim([0 0.79]);

set(gcf,'position',[200 0 160 80]);
set(gca,'position',[0.17 0.33 0.65 0.6]);
clear;