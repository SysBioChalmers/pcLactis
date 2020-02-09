% Data import
load('RcAA_result.mat');

increase_mu = result_rcAA.data(2,:)-result_rcAA.data(1,:);
increase_mu(increase_mu < 0) = 0;
increase_aa = result_rcAA.data(5,:)-result_rcAA.data(3,:);
increase_aa(increase_aa < 0) = 0;
reduced_cost = increase_mu./increase_aa;
scaled_reduced_cost = reduced_cost.*result_rcAA.data(3,:)./result_rcAA.data(1,:);
% scaled_reduced_cost = reduced_cost;
scaled_reduced_cost(isnan(scaled_reduced_cost)) = 0;

aaidlist = result_rcAA.column(1:20);
aaidlist = cellfun(@(x) x(5:end),aaidlist,'UniformOutput',false);

minclr = [255,255,255]/255;
maxclr = [165,15,21]/255;
tmp1 = linspace(minclr(1),maxclr(1),129)';
tmp2 = linspace(minclr(2),maxclr(2),129)';
tmp3 = linspace(minclr(3),maxclr(3),129)';
clrmap = [tmp1 tmp2 tmp3];

figure();
data = [scaled_reduced_cost(1:20);
        scaled_reduced_cost(21:40);
        scaled_reduced_cost(41:60);
        scaled_reduced_cost(61:80)];
h = heatmap(aaidlist,{'refmu0.1','refmu0.3','refmu0.5','refmu0.7'},data,'Colormap',clrmap,'ColorMethod','count','CellLabelColor','none');

set(h,'FontSize',10,'FontName','Helvetica');
% set(gcf,'position',[200 300 100 600]);
% set(gca,'position',[0.2 0.1 0.4 0.75]);


