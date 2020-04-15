load('Sglc_result.mat');
clr = [99,99,99]/255;

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
scatter(x1,y1,10,'o','LineWidth',1,'MarkerEdgeColor',clr,'MarkerEdgeAlpha',0.7);

set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Glucose concentration (uM)','FontSize',7,'FontName','Helvetica');
ylabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');

xlim([0 240]);
ylim([0 0.79]);

set(gcf,'position',[0 0 160 80]);
set(gca,'position',[0.17 0.33 0.55 0.53]);
clear;