% Data import
load('RcAA_result.mat');
[aF, bF, ~] = xlsread('AA_data.xlsx','Flux');
[aC, bC, ~] = xlsread('AA_data.xlsx','Conc');

%% Plot reduced cost analysis
increase_mu = result_rcAA.data(2,:)-result_rcAA.data(1,:);
increase_mu(increase_mu < 0) = 0;
increase_aa = result_rcAA.data(5,:)-result_rcAA.data(3,:);
increase_aa(increase_aa < 0) = 0;
reduced_cost = increase_mu./increase_aa;
scaled_reduced_cost = reduced_cost.*result_rcAA.data(3,:)./result_rcAA.data(1,:);
% scaled_reduced_cost = reduced_cost;
scaled_reduced_cost(isnan(scaled_reduced_cost)) = 0;
% scaled_reduced_cost(scaled_reduced_cost == inf) = 0;

aaidlist = result_rcAA.column(1:20);
aaidlist = cellfun(@(x) x(5:end),aaidlist,'UniformOutput',false);

minclr = [255,255,255]/255;
maxclr = [165,15,21]/255;
tmp1 = linspace(minclr(1),maxclr(1),129)';
tmp2 = linspace(minclr(2),maxclr(2),129)';
tmp3 = linspace(minclr(3),maxclr(3),129)';
clrmap = [tmp1 tmp2 tmp3];

figure('Name','1');
data = [scaled_reduced_cost(1:20);
        scaled_reduced_cost(21:40);
        scaled_reduced_cost(41:60);
        scaled_reduced_cost(61:80)];
h = heatmap(aaidlist,{'refmu0.1','refmu0.3','refmu0.5','refmu0.7'},data,'Colormap',clrmap,'ColorMethod','count','CellLabelColor','none');

set(h,'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 300 350 60]);
set(gca,'position',[0.18 0.3 0.68 0.68]);

%% Plot measured data

flux = zeros(8,length(aaidlist));
conc = zeros(8,length(aaidlist));

for i = 1:length(aaidlist)
    aa_tmp = aaidlist{i};
    if any(contains(bF(1,:),upper(aa_tmp)))
        flux(:,i) = aF(:,contains(bF(1,:),upper(aa_tmp)));
    else
        flux(:,i) = NaN;
    end
    
    if any(contains(bC(1,:),upper(aa_tmp)))
        conc(:,i) = aC(:,contains(bC(1,:),upper(aa_tmp)));
    else
        conc(:,i) = NaN;
    end
end

figure('Name','2');
h = heatmap(conc,'Colormap',clrmap,'ColorMethod','count','CellLabelColor','none');
h.XDisplayLabels = aaidlist;
h.Title = 'Concentration (mM)';
set(h,'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[700 300 350 130]);
set(gca,'position',[0.18 0.3 0.68 0.6]);

figure('Name','3');
h = heatmap(flux,'Colormap',clrmap,'ColorMethod','count','CellLabelColor','none');
h.XDisplayLabels = aaidlist;
h.Title = 'Flux (mmol/gCDW/h)';
set(h,'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[700 500 350 130]);
set(gca,'position',[0.18 0.3 0.68 0.6]);

% figure('Name','4');
% hold on;
% data_ccpa = flux(1:4,:);
% data_wt = flux(5:8,:);
% error_ccpa = std(data_ccpa,1,1);
% error_wt = std(data_wt,1,1);
% 
% b = bar(categorical(aaidlist),[mean(data_ccpa);mean(data_wt)]');
% errorbar([mean(data_ccpa);mean(data_wt)]',[error_ccpa;error_wt]','k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',10);
% 
% boxplot(data_ccpa)




