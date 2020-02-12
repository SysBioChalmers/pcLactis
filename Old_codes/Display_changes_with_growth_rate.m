%% Calculates changes with growth rate. (with sf) 

% Timing: < 10 s

%% Main part
load('Sglc2_fluxes_with_sf.mat');

flux_res = fluxes_simulated_with_sf;

[~, excel_input, ~] = xlsread('Allocation.xlsx');
load('Info_enzyme.mat');
load('Info_mRNA.mat');
load('Info_ribosome.mat');
load('Info_tRNA.mat');

load('pcLactis_Model.mat');
model = pcLactis_Model;
GAM = 36;%ATP coefficient in the new biomass equation.
NGAM = 2; %(mmol/gCDW/h)
f_unmodeled = 0.4; %proportion of unmodeled protein in total protein (g/g)
model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,~] = ChangeUnmodeledProtein(model,f_unmodeled);
unmodeled_weight = 0.46 * f_unmodeled;
clear pcLactis_Model GAM NGAM f_unmodeled;

[~, n] = size(flux_res);

allocation_value = zeros(5,n);
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
    q_glc(1,i) = abs(solME_full(strcmp(model.rxns,'R_M_EX_glc__D_e'),:));
    
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
allocation_pathway = [pathway;'Ribosomal protein';'Inactive enzyme'];
allocation_value = [allocation_value;total_rProtein;total_inactive_enzyme];

allocation_pathway = [allocation_pathway;'Other modeled protein'];
allocation_value = [allocation_value;total_proteome-sum(allocation_value)-unmodeled_weight];

allocation_pathway = [allocation_pathway;'Unmodeled protein'];
allocation_value = [allocation_value;unmodeled_weight*ones(1,n)];

allocation_percentage = allocation_value./total_proteome;


%% Figure
figure('Name','1');
hold on;
b = bar(transpose(mu_list),transpose(allocation_percentage),'stacked');
b(1).FaceColor = [228,26,28]/256;
b(2).FaceColor = [55,126,184]/256;
b(3).FaceColor = [255,255,51]/256;
b(4).FaceColor = [152,78,163]/256;
b(5).FaceColor = [255,127,0]/256;
b(6).FaceColor = [77,175,74]/256;
b(7).FaceColor = [166,86,40]/256;
b(8).FaceColor = [247,129,191]/256;
b(9).FaceColor = [153,153,153]/256;

set(gca,'FontSize',12,'FontName','Helvetica');
legend(allocation_pathway,'FontSize',12,'FontName','Helvetica','location','se');
legend('boxoff');
ylim([0 1]);
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel('Proteome fraction','FontSize',14,'FontName','Helvetica');

set(gcf,'position',[0 0 420 300]);
set(gca,'position',[0.11 0.13 0.5 0.85]);

color = [178,24,43]/255;
% color = [221,28,119]/255;

figure('Name','2');
plot(mu_list,q_glc,'LineWidth',0.75,'Color',color);
ylim([0 28]);
xlim([0 0.75]);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel(['Glucose uptake rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[0 0 270 110]);
set(gca,'position',[0.25 0.28 0.7 0.7]);

figure('Name','3');
plot(mu_list,total_RNA./total_proteome,'LineWidth',0.75,'Color',color);
ylim([0.09 0.24]);
xlim([0 0.75]);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel(['RNA/Protein',char(13,10)','(g/g)'],'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[300 0 270 110]);
set(gca,'position',[0.25 0.28 0.7 0.7]);

figure('Name','4');
plot(mu_list,allocation_value(2,:),'LineWidth',0.75,'Color',color);
ylim([0 0.037]);
xlim([0 0.75]);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel(['Glycolytic enzymes',char(13,10)','(g/gCDW)'],'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[600 0 270 110]);
set(gca,'position',[0.25 0.28 0.7 0.7]);

figure('Name','5');
plot(mu_list,total_inactive_enzyme,'LineWidth',0.75,'Color',color);
ylim([0 0.25]);
xlim([0 0.75]);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel(['Inactive enzymes',char(13,10)','(g/gCDW)'],'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[900 0 270 110]);
set(gca,'position',[0.25 0.28 0.7 0.7]);