
load('Sglc_fluxes');
% load('Sglc_fluxes_test.mat');
flux_res = fluxes_simulated_without_sf;

load('pcLactis_Model.mat');
model = pcLactis_Model;
GAM = 36; %ATP coefficient in the new biomass equation.
NGAM = 2; %(mmol/gCDW/h)
f_unmodeled = 0.4; %proportion of unmodeled protein in total protein (g/g)
model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,~] = ChangeUnmodeledProtein(model,f_unmodeled);
unmodeled_weight = 0.46 * f_unmodeled;
clear pcLactis_Model GAM NGAM f_unmodeled;

%% Metabolic shift
% Experimental data
[~, ~, exp_raw] = xlsread('Exchange_reaction_setting.xlsx','Exchange_rates');
exp_mu = cell2mat(exp_raw(2,14:18));

exp_arg = 1*cell2mat(exp_raw(9,14:18));
exp_orn = cell2mat(exp_raw(8,14:18));
exp_nh4 = cell2mat(exp_raw(10,14:18));
exp_citr = cell2mat(exp_raw(11,14:18));
sd_orn = cell2mat(exp_raw(8,19:23));
sd_arg = cell2mat(exp_raw(9,19:23));
sd_nh4 = cell2mat(exp_raw(10,19:23));
sd_citr = cell2mat(exp_raw(11,19:23));

mu = flux_res(strcmp(model.rxns,'R_biomass_dilution'),:);
orn = flux_res(strcmp(model.rxns,'R_M_EX_orn_e'),:);
nh4 = flux_res(strcmp(model.rxns,'R_M_EX_nh4_e'),:);
arg = flux_res(strcmp(model.rxns,'R_M_EX_arg__L_e'),:);
citr = flux_res(strcmp(model.rxns,'R_M_EX_citr__L_e'),:);


color_arg = [231,41,138]/255;
color_orn = [27,158,119]/255;
color_nh4 = [217,95,2]/255;
color_citr = [117,112,179]/255;

figure('Name','1');

hold on;
box on;

scatter(exp_mu,exp_arg,10,'o','LineWidth',1,'MarkerEdgeColor',color_arg);
scatter(exp_mu,exp_orn,10,'o','LineWidth',1,'MarkerEdgeColor',color_orn);
scatter(exp_mu,exp_nh4,10,'o','LineWidth',1,'MarkerEdgeColor',color_nh4);
scatter(exp_mu,exp_citr,10,'o','LineWidth',1,'MarkerEdgeColor',color_citr);

plot(mu,arg,'-','Marker','.','MarkerSize',5,'LineWidth',0.75,'Color',color_arg);
plot(mu,orn,'-','Marker','.','MarkerSize',5,'LineWidth',0.75,'Color',color_orn);
plot(mu,nh4,'-','Marker','.','MarkerSize',5,'LineWidth',0.75,'Color',color_nh4);
plot(mu,citr,'-','Marker','.','MarkerSize',5,'LineWidth',0.75,'Color',color_citr);

xlim([0.1 0.7]);
% ylim([0 1.5]);
set(gca, 'XColor','k');
set(gca, 'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica','Color','k');
ylabel('Exchange flux (mmol/gCDW/h)','FontSize',7,'FontName','Helvetica','Color','k');
legend({'Arginine','Ornithine','Ammonium','﻿Citrulline',...
    'Arginine','Ornithine','Ammonium','﻿Citrulline'},'FontSize',6,'FontName','Helvetica','location','se');
title('Arginine catabolism','FontSize',7,'FontName','Helvetica','Color','k');

set(gcf,'position',[200 400 300 200]);
set(gca,'position',[0.15 0.2 0.4 0.7]);


%% Protein constraints
[~, excel_input, ~] = xlsread('Allocation.xlsx');
load('Info_enzyme.mat');
load('Info_mRNA.mat');
load('Info_ribosome.mat');
load('Info_tRNA.mat');

[~, n] = size(flux_res);

allocation_value = zeros(5,n);
total_proteome = zeros(1,n);
total_RNA = zeros(1,n);
total_rProtein = zeros(1,n);
total_inactive_enzyme = zeros(1,n);
mu_list = zeros(1,n);
total_Ribo = zeros(1,n);
q_glc = zeros(1,n);
total_rRNA = zeros(1,n);
total_tRNA = zeros(1,n);
total_mRNA = zeros(1,n);

for i = 1:n
    
    disp([num2str(i) '/' num2str(n)]);
    
    solME_full = flux_res(:,i);
    
    mu = solME_full(strcmp(model.rxns,'R_biomass_dilution'));
    mu_list(1,i) = mu;
    q_glc(1,i) = abs(solME_full(strcmp(model.rxns,'R_M_EX_glc__D_e'),:));
    
    [pathway, absvalue] = CalculateProteinAllocation(model,solME_full,Info_enzyme,excel_input);
    allocation_value(:,i) = absvalue;
    
    [Protein,RNA,rProtein,rRNA,tRNA,mRNA] = CalculateProteinAndRNA(model,solME_full,...
                                                             Info_enzyme,...
                                                             Info_ribosome,...
                                                             Info_tRNA,...
                                                             Info_mRNA);
    total_proteome(1,i) = Protein;
    total_RNA(1,i) = RNA;
    total_rProtein(1,i) = rProtein;
    total_Ribo(1,i) = rProtein + rRNA;
    
    factor_k = 1;
    
    [~,~,f_enzyme_inact,~] = CheckInactiveEnzyme(model,solME_full,factor_k);
    total_inactive_enzyme(1,i) = f_enzyme_inact;
    total_rRNA(1,i) = rRNA;
    total_tRNA(1,i) = tRNA;
    total_mRNA(1,i) = mRNA;
end
allocation_pathway = [pathway;'Ribosomal protein';'Inactive enzyme'];
allocation_value = [allocation_value;total_rProtein;total_inactive_enzyme];

allocation_pathway = [allocation_pathway;'Other modelled protein'];
allocation_value = [allocation_value;total_proteome-sum(allocation_value)-unmodeled_weight];

allocation_pathway = [allocation_pathway;'Unmodelled protein'];
allocation_value = [allocation_value;unmodeled_weight*ones(1,n)];

allocation_percentage = allocation_value./total_proteome;

% Figures
figure('Name','2');
hold on;
b = bar(transpose(allocation_percentage),'stacked');
b(1).FaceColor = [228,26,28]/256;
b(2).FaceColor = [55,126,184]/256;
b(3).FaceColor = [255,255,51]/256;
b(4).FaceColor = [152,78,163]/256;
b(5).FaceColor = [255,127,0]/256;
b(6).FaceColor = [77,175,74]/256;
b(7).FaceColor = [166,86,40]/256;
b(8).FaceColor = [247,129,191]/256;
b(9).FaceColor = [153,153,153]/256;


legend(allocation_pathway,'FontSize',6,'FontName','Helvetica','location','se');
legend('boxoff');
ylim([0 1]);
xticks(1:1:length(mu_list));
xticklabels(cellfun(@(x) num2str(x),num2cell(mu_list),'UniformOutput',false));
xtickangle(90);
xlabel('Growth rate (/h)','FontSize',6,'FontName','Helvetica');
ylabel('Proteome fraction','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
set(gcf,'position',[200 100 420 200]);
set(gca,'position',[0.1 0.2 0.5 0.7]);




