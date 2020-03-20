%% Plot growth rate and saturation factor

load('Egsf2_result.mat');
clr = [99,99,99]/255;
x = global_saturation_factor_list(:,1);
y = global_saturation_factor_list(:,2);
% x = x(y ~= 1);
% y = y(y ~= 1);
x = x(~isnan(y));
y = y(~isnan(y));
p = polyfit(x(1:end-1),y(1:end-1),1);

hold on;
box on;
scatter(x,y,20,'x','LineWidth',0.1,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.7);

x1 = linspace(x(1),x(end-1));  
y1 = polyval(p,x1); 

plot(x1,y1,':','LineWidth',1,'Color',clr);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
ylabel('Global saturation','FontSize',7,'FontName','Helvetica');

equation = ['y = ',num2str(round(p(1),3)),'x',num2str(round(p(2),3)),];
text(0.25,0.25,equation,'FontSize',7,'FontName','Helvetica','Color','k');

rsq = corr(y(1:end-1),polyval(p,x(1:end-1)))^2;
R2 = ['R^2 = ',num2str(round(rsq,3))];
text(0.12,0.8,R2,'FontSize',7,'FontName','Helvetica','Color','k');

ylim([0 1]);
xlim([0.1 0.6]);
set(gca,'xtick',0.1:0.1:0.6);

set(gcf,'position',[0 0 160 80]);
set(gca,'position',[0.17 0.33 0.55 0.53]);


% load('Egsf2_result.mat');
% clr = [99,99,99]/255;
% x = global_saturation_factor_list(:,1);
% y = global_saturation_factor_list(:,2);
% % x = x(y ~= 1);
% % y = y(y ~= 1);
% x = x(~isnan(y));
% y = y(~isnan(y));
% k = x\y;
% 
% hold on;
% box on;
% scatter(x,y,20,'x','LineWidth',0.1,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.7);
% plot([0;0.6],[0;0.6*k],':','LineWidth',1,'Color',clr);
% set(gca,'FontSize',6,'FontName','Helvetica');
% xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
% ylabel('Global saturation factor','FontSize',7,'FontName','Helvetica');
% 
% equation = ['y = ',num2str(round(k,3)),' x'];
% text(0.05,0.8,equation,'FontSize',7,'FontName','Helvetica','Color','k');
% 
% yfit =  k * x;
% yresid = y - yfit;
% SSresid = sum(yresid.^2);
% SStotal = (length(y)-1) * var(y);
% rsq = 1 - SSresid/SStotal;
% 
% R2 = ['R^2 = ',num2str(round(rsq,3))];
% text(0.05,0.6,R2,'FontSize',7,'FontName','Helvetica','Color','k');
% 
% ylim([0 1]);
% 
% set(gcf,'position',[0 0 160 80]);
% set(gca,'position',[0.17 0.33 0.55 0.53]);

