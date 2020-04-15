%% Figure
load('Egt_result.mat');
figure('Name','1');
hold on;
box on;
x = res(:,9);
y = res(:,7);
plot(x,y,'-o','LineWidth',1,'MarkerSize',1);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Fraction of glucose transporter','FontSize',8,'FontName','Helvetica');
ylabel('Growth rate (/h)','FontSize',8,'FontName','Helvetica');

% xlim([0 0.2]);
ylim([0 0.8]);


set(gcf,'position',[200 0 160 80]);
set(gca,'position',[0.17 0.33 0.75 0.53]);

clear x y;
