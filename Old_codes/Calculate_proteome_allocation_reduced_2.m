%% This function calculates proteome allocation for exp data.

% After changing saturation factor.

load('Egsf1_fluxes.mat');

flux_res = fluxes_global_saturation_factor_changed;

[~, excel_input, ~] = xlsread('Allocation.xlsx');
excel_input = excel_input(:,2);
load('Info_enzyme.mat');
load('Info_mRNA.mat');
load('Info_ribosome.mat');
load('Info_tRNA.mat');

load('pcLactis_Model.mat');
model = pcLactis_Model;
GAM = 42;%ATP coefficient in the new biomass equation.
NGAM = 2.5; %(mmol/gCDW/h)
f_unmodeled = 0.42; %proportion of unmodeled protein in total protein (g/g)
model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,~] = ChangeUnmodeledProtein(model,f_unmodeled);
unmodeled_weight = 0.46 * f_unmodeled;
clear pcLactis_Model GAM NGAM f_unmodeled;

[~, n] = size(flux_res);

allocation_value = zeros(1,n);
total_proteome = zeros(1,n);
total_RNA = zeros(1,n);
total_rProtein = zeros(1,n);
total_inactive_enzyme = zeros(1,n);
mu_list = zeros(1,n);
total_Ribo = zeros(1,n);
q_glc = zeros(1,n);

for i = 1:n
    
    disp([num2str(i) '/' num2str(n)]);
    
    solME_full = flux_res(:,i);
    
    mu = solME_full(strcmp(model.rxns,'R_biomass_dilution'));
    mu_list(1,i) = mu;
    q_glc(1,i) = abs(solME_full(strcmp(model.rxns,'R_M_EX_glc_LPAREN_e_RPAREN_'),:));
    
    [pathway, absvalue] = CalculateProteinAllocation(model,solME_full,Info_enzyme,excel_input);
    allocation_value(:,i) = absvalue;
    
    [Protein,RNA,rProtein,rRNA,~,~] = CalculateProteinAndRNA(model,solME_full,...
                                                             Info_enzyme,...
                                                             Info_ribosome,...
                                                             Info_tRNA,...
                                                             Info_mRNA);
    total_proteome(1,i) = Protein;
    total_RNA(1,i) = RNA;
    total_rProtein(1,i) = rProtein;
    total_Ribo(1,i) = rProtein + rRNA;
    [~,~,f_enzyme_inact] = CheckInactiveEnzyme(model,solME_full);
    total_inactive_enzyme(1,i) = f_enzyme_inact;
    
end
allocation_pathway = [pathway;'Ribosomal proteins';'Inactive enzymes'];
allocation_value = [allocation_value;total_rProtein;total_inactive_enzyme];

allocation_pathway = [allocation_pathway;'Other modeled proteins'];
allocation_value = [allocation_value;total_proteome-sum(allocation_value)-unmodeled_weight];

allocation_pathway = [allocation_pathway;'Unmodeled proteins'];
allocation_value = [allocation_value;unmodeled_weight*ones(1,n)];

allocation_percentage = allocation_value./total_proteome;


%% Figure

perc_015 = allocation_percentage(:,1:3);
perc_030 = allocation_percentage(:,4:6);
perc_045 = allocation_percentage(:,7);
perc_050 = allocation_percentage(:,8:9);
perc_060 = allocation_percentage(:,10);
mean_015 = mean(perc_015,2);
mean_030 = mean(perc_030,2);
mean_045 = mean(perc_045,2);
mean_050 = mean(perc_050,2);
mean_060 = mean(perc_060,2);
error_015 = std(perc_015,1,2);
error_030 = std(perc_030,1,2);
error_045 = std(perc_045,1,2);
error_050 = std(perc_050,1,2);
error_060 = std(perc_060,1,2);

res_mean = [mean_015 mean_030 mean_045 mean_050 mean_060];
res_error = [error_015 error_030 error_045 error_050 error_060];
x = categorical({'0.15','0.30','0.45','0.50','0.60',});

figure('Name','1');
hold on;
b = bar(x,res_mean','stacked');

b(4).FaceColor = [253,191,111]/255;
b(1).FaceColor = [178,223,138]/255;
b(2).FaceColor = [251,154,153]/255;
b(3).FaceColor = [166,206,227]/255;
b(5).FaceColor = [178,178,178]/255;
errorbar(cumsum(res_mean)',res_error','k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);

set(gca,'FontSize',12,'FontName','Helvetica');
legend(allocation_pathway,'FontSize',12,'FontName','Helvetica','location','se');
legend('boxoff');
ylim([0 1]);
ylabel('Proteome fraction','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');

set(gcf,'position',[0 0 250 300]);
set(gca,'position',[0.2 0.13 0.7 0.51]);