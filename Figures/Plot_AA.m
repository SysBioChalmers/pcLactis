% Data import
load('RcAA_result.mat');
load('Sglc_fluxes');
% load('Sglc_fluxes_test.mat');
flux_res = fluxes_simulated_without_sf;

load('pcLactis_Model.mat');
model = pcLactis_Model;

%% Plot reduced cost analysis

increase_mu = result_rcAA.data(2,:)-result_rcAA.data(1,:);
increase_mu(increase_mu < 0) = 0;
increase_aa = result_rcAA.data(5,:)-result_rcAA.data(3,:);
increase_aa(increase_aa < 0) = 0;
reduced_cost = increase_mu./increase_aa;
scaled_reduced_cost = reduced_cost.*result_rcAA.data(3,:)./result_rcAA.data(1,:);

tf_aa = result_rcAA.data(5,:)-result_rcAA.data(4,:);
tf_aa(abs(tf_aa)<1e-12) = 0;
scaled_reduced_cost(tf_aa ~= 0) = 0;
scaled_reduced_cost(isnan(scaled_reduced_cost)) = 0;
% scaled_reduced_cost(scaled_reduced_cost == inf) = 0;

aaidlist = result_rcAA.column(1:20);
aaidlist = cellfun(@(x) x(5:end),aaidlist,'UniformOutput',false);

mulist = result_rcAA.column;
mulist = cellfun(@(x) x(1:strfind(x,'_')-1),mulist,'UniformOutput',false);
mulist = unique(mulist);

minclr = [255,255,255]/255;
maxclr = [197,27,138]/255;
tmp1 = linspace(minclr(1),maxclr(1),129)';
tmp2 = linspace(minclr(2),maxclr(2),129)';
tmp3 = linspace(minclr(3),maxclr(3),129)';
clrmap = [tmp1 tmp2 tmp3];

figure('Name','1');
data = [scaled_reduced_cost(1:20);
        scaled_reduced_cost(21:40);
        scaled_reduced_cost(41:60);
        scaled_reduced_cost(61:80)];
h = heatmap(mulist,aaidlist,data','Colormap',clrmap,'ColorMethod','count','CellLabelColor','none');
h.XLabel = 'Growth rate (/h)';
title('Scaled reduced cost');
set(h,'FontSize',6.5,'FontName','Helvetica');
h.FontColor = 'k';
set(gcf,'position',[400 300 160 320]);
set(gca,'position',[0.13 0.13 0.35 0.6]);


%% Simulated AA
[~, ~, aa_raw] = xlsread('Exchange_reaction_setting.xlsx','AA_factors');
Exchange_AAs = aa_raw(2:end,1);
LBfactor_AAs = cell2mat(aa_raw(2:end,2));
clear aa_raw;

[~, ~, aa_data] = xlsread('Exchange_reaction_setting.xlsx','chemostat_data_2');
aa_list = aa_data(1,2:end);
aa_value = cell2mat(aa_data(2:end,2:end));
mu_exp = cell2mat(aa_data(2:end,1));
clear aa_data;

mu2 = flux_res(strcmp(model.rxns,'R_biomass_dilution'),:);

% clr = [69,117,180]/255;
clr = maxclr;
figure('Name','2');
for i = 1:length(Exchange_AAs)
    rxnid = Exchange_AAs{i};
    aaid = rxnid(8:10);
    flux_tmp = -flux_res(strcmp(model.rxns,rxnid),:);
    flux_tRNA = CalculateAAtRNAFlux(model,flux_res,aaid);
    m = -LBfactor_AAs(i);
    
    subplot(4,5,i);
    hold on;
    box on;
    plot(mu2,flux_tmp,'-','LineWidth',1,'Color',clr);
    plot(mu2,flux_tRNA,'-','LineWidth',1,'Color','k');
    x = [0,1]; y = [0,m];
    plot(x,y,':','LineWidth',1.5,'Color','k');
    
    if contains(aaid,aa_list)
        aa_exp = -aa_value(:,contains(aa_list,aaid));
        scatter(mu_exp,aa_exp,15,'o','filled',...
            'LineWidth',1,'MarkerEdgeColor',[1 1 1],...
            'MarkerFaceColor',clr,...
            'MarkerEdgeAlpha',0,...
            'MarkerFaceAlpha',0.5);
    end
    
    xlim([0.1 0.7]);
    set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
%     ylabel('Uptake flux','FontSize',14,'FontName','Helvetica');
%     xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
    title(aaid,'FontSize',7,'FontName','Helvetica','Color','k');
end

set(gcf,'position',[0 0 350 280]);
