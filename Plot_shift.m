load('Sglc2_fluxes_with_sf.mat');
flux_res = fluxes_simulated_with_sf;

load('pcLactis_Model.mat');
model = pcLactis_Model;
GAM = 38;%ATP coefficient in the new biomass equation.
NGAM = 2; %(mmol/gCDW/h)
f_unmodeled = 0.45; %proportion of unmodeled protein in total protein (g/g)
model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,~] = ChangeUnmodeledProtein(model,f_unmodeled);
unmodeled_weight = 0.46 * f_unmodeled;
clear pcLactis_Model GAM NGAM f_unmodeled;

%% Metabolic shift
% Experimental data
[~, ~, exp_raw] = xlsread('Exchange_reaction_setting.xlsx','Exchange_rates');
exp_mu = cell2mat(exp_raw(2,14:18));

exp_glc = -1*cell2mat(exp_raw(3,14:18));
exp_ac = cell2mat(exp_raw(4,14:18));
exp_eth = cell2mat(exp_raw(5,14:18));
exp_form = cell2mat(exp_raw(6,14:18));
exp_lac = cell2mat(exp_raw(7,14:18));
sd_glc = cell2mat(exp_raw(3,19:23));
sd_ac = cell2mat(exp_raw(4,19:23));
sd_eth = cell2mat(exp_raw(5,19:23));
sd_form = cell2mat(exp_raw(6,19:23));
sd_lac = cell2mat(exp_raw(7,19:23));

exp_arg = -1*cell2mat(exp_raw(9,14:18));
exp_orn = cell2mat(exp_raw(8,14:18));
exp_nh4 = cell2mat(exp_raw(10,14:18));
exp_citr = cell2mat(exp_raw(11,14:18));
sd_orn = cell2mat(exp_raw(8,19:23));
sd_arg = cell2mat(exp_raw(9,19:23));
sd_nh4 = cell2mat(exp_raw(10,19:23));
sd_citr = cell2mat(exp_raw(11,19:23));

mu2 = flux_res(strcmp(model.rxns,'R_biomass_dilution'),:);

glc2 = -flux_res(strcmp(model.rxns,'R_M_EX_glc__D_e'),:);
ac2 = flux_res(strcmp(model.rxns,'R_M_EX_ac_e'),:);
eth2 = flux_res(strcmp(model.rxns,'R_M_EX_etoh_e'),:);
form2 = flux_res(strcmp(model.rxns,'R_M_EX_for_e'),:);
lac2 = flux_res(strcmp(model.rxns,'R_M_EX_lac__L_e'),:);

orn2 = flux_res(strcmp(model.rxns,'R_M_EX_orn_e'),:);
citr2 = flux_res(strcmp(model.rxns,'R_M_EX_citr__L_e'),:);
nh42 = flux_res(strcmp(model.rxns,'R_M_EX_nh4_e'),:);
arg2 = -flux_res(strcmp(model.rxns,'R_M_EX_arg__L_e'),:);


% idx = [1:5,8:13];
% mu2 = mu2(idx);
% glc2 = glc2(idx);
% ac2 = ac2(idx);
% eth2 = eth2(idx);
% form2 = form2(idx);
% lac2 = lac2(idx);

% Figures
color_glc = [228,26,28]/255;
color_ac = [55,126,184]/255;
color_eth = [255,127,0]/255;
color_form = [77,175,74]/255;
color_lac = [152,78,163]/255;

color_arg = [247,129,191]/255;
color_orn = [166,86,40]/255;
color_nh4 = [153,153,153]/255;
color_citr = [0,0,0]/255;

figure('Name','1');
subplot(2,5,1);
hold on;
box on;
plot(mu2,glc2,'-','LineWidth',0.75,'Color',color_glc);
plot(exp_mu,exp_glc,'-.','LineWidth',0.75,'Color',color_glc);
scatter(mu2,glc2,10,'o','filled','LineWidth',0.75,'MarkerEdgeColor',color_glc,'MarkerFaceColor',color_glc);
scatter(exp_mu,exp_glc,10,'o','LineWidth',0.75,'MarkerEdgeColor',color_glc);
x = exp_mu; y = exp_glc; yu = exp_glc + sd_glc; yl = exp_glc - sd_glc;
plot(x,y,'-.','LineWidth',0.75,'Color',color_glc);
fill([x fliplr(x)],[yu fliplr(yl)],color_glc,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Flux','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('Glucose','FontSize',14,'FontName','Helvetica');

subplot(2,5,2);
hold on;
box on;
plot(mu2,ac2,'-','LineWidth',0.75,'Color',color_ac);
plot(exp_mu,exp_ac,'-.','LineWidth',0.75,'Color',color_ac);
scatter(mu2,ac2,10,'o','filled','LineWidth',0.75,'MarkerEdgeColor',color_ac,'MarkerFaceColor',color_ac);
scatter(exp_mu,exp_ac,10,'o','LineWidth',0.75,'MarkerEdgeColor',color_ac);
x = exp_mu; y = exp_ac; yu = exp_ac + sd_ac; yl = exp_ac - sd_ac;
plot(x,y,'-.','LineWidth',0.75,'Color',color_ac);
fill([x fliplr(x)],[yu fliplr(yl)],color_ac,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Flux','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('Acetate','FontSize',14,'FontName','Helvetica');

subplot(2,5,3);
hold on;
box on;
plot(mu2,eth2,'-','LineWidth',0.75,'Color',color_eth);
plot(exp_mu,exp_eth,'-.','LineWidth',0.75,'Color',color_eth);
scatter(mu2,eth2,10,'o','filled','LineWidth',0.75,'MarkerEdgeColor',color_eth,'MarkerFaceColor',color_eth);
scatter(exp_mu,exp_eth,10,'o','LineWidth',0.75,'MarkerEdgeColor',color_eth);
x = exp_mu; y = exp_eth; yu = exp_eth + sd_eth; yl = exp_eth - sd_eth;
plot(x,y,'-.','LineWidth',0.75,'Color',color_eth);
fill([x fliplr(x)],[yu fliplr(yl)],color_eth,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Flux','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('Ethanol','FontSize',14,'FontName','Helvetica');

subplot(2,5,4);
hold on;
box on;
plot(mu2,form2,'-','LineWidth',0.75,'Color',color_form);
plot(exp_mu,exp_form,'-.','LineWidth',0.75,'Color',color_form);
scatter(mu2,form2,10,'o','filled','LineWidth',0.75,'MarkerEdgeColor',color_form,'MarkerFaceColor',color_form);
scatter(exp_mu,exp_form,10,'o','LineWidth',0.75,'MarkerEdgeColor',color_form);
x = exp_mu; y = exp_form; yu = exp_form + sd_form; yl = exp_form - sd_form;
plot(x,y,'-.','LineWidth',0.75,'Color',color_form);
fill([x fliplr(x)],[yu fliplr(yl)],color_form,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Flux','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('Formate','FontSize',14,'FontName','Helvetica');

subplot(2,5,5);
hold on;
box on;
plot(mu2,lac2,'-','LineWidth',0.75,'Color',color_lac);
plot(exp_mu,exp_lac,'-.','LineWidth',0.75,'Color',color_lac);
scatter(mu2,lac2,10,'o','filled','LineWidth',0.75,'MarkerEdgeColor',color_lac,'MarkerFaceColor',color_lac);
scatter(exp_mu,exp_lac,10,'o','LineWidth',0.75,'MarkerEdgeColor',color_lac);
x = exp_mu; y = exp_lac; yu = exp_lac + sd_lac; yl = exp_lac - sd_lac;
plot(x,y,'-.','LineWidth',0.75,'Color',color_lac);
fill([x fliplr(x)],[yu fliplr(yl)],color_lac,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Flux','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('Lactate','FontSize',14,'FontName','Helvetica');

subplot(2,5,6);
hold on;
box on;
plot(mu2,arg2,'-','LineWidth',0.75,'Color',color_arg);
plot(exp_mu,exp_arg,'-.','LineWidth',0.75,'Color',color_arg);
scatter(mu2,arg2,10,'o','filled','LineWidth',0.75,'MarkerEdgeColor',color_arg,'MarkerFaceColor',color_arg);
scatter(exp_mu,exp_arg,10,'o','LineWidth',0.75,'MarkerEdgeColor',color_arg);
x = exp_mu; y = exp_arg; yu = exp_arg + sd_arg; yl = exp_arg - sd_arg;
plot(x,y,'-.','LineWidth',0.75,'Color',color_arg);
fill([x fliplr(x)],[yu fliplr(yl)],color_arg,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Flux','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('Arginine','FontSize',14,'FontName','Helvetica');

subplot(2,5,7);
hold on;
box on;
plot(mu2,orn2,'-','LineWidth',0.75,'Color',color_orn);
plot(exp_mu,exp_orn,'-.','LineWidth',0.75,'Color',color_orn);
scatter(mu2,orn2,10,'o','filled','LineWidth',0.75,'MarkerEdgeColor',color_orn,'MarkerFaceColor',color_orn);
scatter(exp_mu,exp_orn,10,'o','LineWidth',0.75,'MarkerEdgeColor',color_orn);
x = exp_mu; y = exp_orn; yu = exp_orn + sd_orn; yl = exp_orn - sd_orn;
plot(x,y,'-.','LineWidth',0.75,'Color',color_orn);
fill([x fliplr(x)],[yu fliplr(yl)],color_orn,'linestyle','none','FaceAlpha',0.3);
xlim([0.1 0.7]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Flux','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
title('Ornithine','FontSize',14,'FontName','Helvetica');

% subplot(2,5,8);
% hold on;
% box on;
% plot(mu2,nh42,'-','LineWidth',0.75,'Color',color_nh4);
% plot(exp_mu,exp_nh4,'-.','LineWidth',0.75,'Color',color_nh4);
% xlim([0.1 0.7]);
% set(gca,'FontSize',12,'FontName','Helvetica');
% ylabel('Flux','FontSize',14,'FontName','Helvetica');
% xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
% title('Ammonium','FontSize',14,'FontName','Helvetica');
% 
% subplot(2,5,9);
% hold on;
% box on;
% plot(mu2,citr2,'-','LineWidth',0.75,'Color',color_citr);
% plot(exp_mu,exp_citr,'-.','LineWidth',0.75,'Color',color_citr);
% xlim([0.1 0.7]);
% set(gca,'FontSize',12,'FontName','Helvetica');
% ylabel('Flux','FontSize',14,'FontName','Helvetica');
% xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
% title('Citrulline','FontSize',14,'FontName','Helvetica');

%% AA constraints
[~, ~, aa_raw] = xlsread('Exchange_reaction_setting.xlsx','AA_factors');
Exchange_AAs = aa_raw(2:end,1);
LBfactor_AAs = cell2mat(aa_raw(2:end,2));
clear aa_raw;

[~, ~, aa_data] = xlsread('Exchange_reaction_setting.xlsx','chemostat_data_2');
aa_list = aa_data(1,2:end);
aa_value = cell2mat(aa_data(2:end,2:end));
mu_exp = cell2mat(aa_data(2:end,1));
clear aa_data;

figure('Name','2');
for i = 1:length(Exchange_AAs)
    rxnid = Exchange_AAs{i};
    aaid = rxnid(8:10);
    flux_tmp = -flux_res(strcmp(model.rxns,rxnid),:);
    m = -LBfactor_AAs(i);
    
    subplot(4,5,i);
    hold on;
    box on;
    plot(mu2,flux_tmp,'-','LineWidth',0.75,'Color',[247,104,161]/255);
    x = [0,1]; y = [0,m];
    plot(x,y,':','LineWidth',1.5,'Color','k');
    
    if contains(aaid,aa_list)
        aa_exp = -aa_value(:,contains(aa_list,aaid));
        scatter(mu_exp,aa_exp,30,'o','LineWidth',0.75,'MarkerEdgeColor',[197,27,138]/255);
    end
    
    xlim([0.1 0.7]);
    set(gca,'FontSize',12,'FontName','Helvetica');
    ylabel('Uptake flux','FontSize',14,'FontName','Helvetica');
    xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
    title(aaid,'FontSize',14,'FontName','Helvetica');
end

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

% Figures
figure('Name','3');
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

figure('Name','4');
plot(mu_list,q_glc,'LineWidth',0.75,'Color',color);
ylim([0 28]);
xlim([0 0.75]);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel(['Glucose uptake rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[0 0 270 110]);
set(gca,'position',[0.25 0.28 0.7 0.7]);

figure('Name','5');
plot(mu_list,total_RNA./total_proteome,'LineWidth',0.75,'Color',color);
% ylim([0.1 0.24]);
xlim([0 0.75]);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel(['RNA/Protein',char(13,10)','(g/g)'],'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[300 0 270 110]);
set(gca,'position',[0.25 0.28 0.7 0.7]);

figure('Name','6');
plot(mu_list,allocation_value(2,:),'LineWidth',0.75,'Color',color);
% ylim([0 0.037]);
xlim([0 0.75]);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel(['Glycolytic enzymes',char(13,10)','(g/gCDW)'],'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[600 0 270 110]);
set(gca,'position',[0.25 0.28 0.7 0.7]);

figure('Name','7');
plot(mu_list,total_inactive_enzyme,'LineWidth',0.75,'Color',color);
ylim([0 0.012]);
xlim([0 0.75]);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel(['Inactive enzymes',char(13,10)','(g/gCDW)'],'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[900 100 270 170]);
set(gca,'position',[0.25 0.28 0.7 0.5]);

figure('Name','8');
hold on;
box on;
x = [0,1]; y = [0.0046,0.0046];
plot(x,y,':','LineWidth',1.5,'Color','k');
plot(mu_list,allocation_value(1,:),'LineWidth',0.75,'Color',color);
ylim([0.004 0.005]);
xlim([0 0.75]);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel(['Glucose transporter',char(13,10)','(g/gCDW)'],'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[900 200 270 170]);
set(gca,'position',[0.25 0.28 0.7 0.5]);


