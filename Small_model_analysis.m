% Here will calculate protein efficiency and ATP yield for mixed acid
% pathway, lactic acid pathway and arginine catabolism.

kcat_glc = 180;%kcat value of glucose transporter % /s

% Load model
load('M_Model.mat');
model = M_Model;
clear M_Model;
model.lb(contains(model.rxns,'R_M_EX_')) = 0;
model.ub(contains(model.rxns,'R_M_EX_')) = 1000;
model = changeRxnBounds(model,'R_M_ATPM',0,'l');
model = changeRxnBounds(model,'R_M_ATPM',1000,'u');
model.lb(contains(model.rxns,'R_M_biomass_LLA')) = 0;
model.ub(contains(model.rxns,'R_M_biomass_LLA')) = 0;
model = changeRxnBounds(model,'R_M_GLCt2',0,'b');

%% Glycolysis + mixed acid fermentation
model_1 = changeObjective(model,'R_M_ATPM',1);
model_1 = changeRxnBounds(model_1,'R_M_EX_glc__D_e',-1,'b');
model_1 = changeRxnBounds(model_1,'R_M_EX_h2o_e',-1000,'l');
sol_1 = optimizeCbModel(model_1,'max','one');
rxnlist_1 = model_1.rxns(sol_1.x ~= 0);
fluxlist_1 = sol_1.x(sol_1.x ~= 0);
yatp_1 = sol_1.x(ismember(model_1.rxns,'R_M_ATPM')) / -sol_1.x(ismember(model_1.rxns,'R_M_EX_glc__D_e'));
clear model_1 sol_1;

%% Glycolysis + lactic acid fermentation
model_2 = changeObjective(model,'R_M_EX_lac__L_e',1);
model_2 = changeRxnBounds(model_2,'R_M_EX_glc__D_e',-1,'b');
model_2 = changeRxnBounds(model_2,'R_M_ATPM',2,'b');
sol_2 = optimizeCbModel(model_2,'max','one');
rxnlist_2 = model_2.rxns(sol_2.x ~= 0);
fluxlist_2 = sol_2.x(sol_2.x ~= 0);
yatp_2 = sol_2.x(ismember(model_2.rxns,'R_M_ATPM')) / -sol_2.x(ismember(model_2.rxns,'R_M_EX_glc__D_e'));
clear model_2 sol_2;

%% Arginine catabolism
model_3 = changeObjective(model,'R_M_ATPM',1);
model_3 = changeRxnBounds(model_3,'R_M_EX_arg__L_e',-1,'b');
model_3 = changeRxnBounds(model_3,'R_M_EX_h2o_e',-1000,'l');
model_3 = changeRxnBounds(model_3,'R_M_ARGt2r',0,'b');
model_3 = changeRxnBounds(model_3,'R_M_ARGabc',0,'b');
model_3 = changeRxnBounds(model_3,'R_M_EX_h_e',-1000,'l');
sol_3 = optimizeCbModel(model_3,'max','one');
rxnlist_3 = model_3.rxns(sol_3.x ~= 0);
fluxlist_3 = sol_3.x(sol_3.x ~= 0);
yatp_3 = sol_3.x(ismember(model_3.rxns,'R_M_ATPM')) / -sol_3.x(ismember(model_3.rxns,'R_M_EX_arg__L_e'));
clear model_3 sol_3;

%% Minimal protein cost for each rxn and pathway
load('Info_enzyme.mat');
[num, txt, ~] = xlsread('k_parameter.xlsx','M');
kcatlist_tmp = num;
enzymelist_tmp = txt(2:end,1);
kcatlist_tmp(kcatlist_tmp < quantile(kcatlist_tmp,0.1,1)) = quantile(kcatlist_tmp,0.1,1);
% remove one of ADH isozymes llmg_0955
enzymelist = enzymelist_tmp(~ismember(enzymelist_tmp,'M_ALCD2x_1_rvs_Enzyme_c'));
kcatlist = kcatlist_tmp(~ismember(enzymelist_tmp,'M_ALCD2x_1_rvs_Enzyme_c'));
clear num txt enzymelist_tmp kcatlist_tmp;

% Glycolysis + mixed acid fermentation
[protcost_1,totprotcost_1] = CalculateProteinCost(rxnlist_1,fluxlist_1,Info_enzyme,kcat_glc,enzymelist,kcatlist);
protein_efficiency_1 = yatp_1/totprotcost_1; %mmolATP/gProtein/h

% Glycolysis + lactic acid fermentation
[protcost_2,totprotcost_2] = CalculateProteinCost(rxnlist_2,fluxlist_2,Info_enzyme,kcat_glc,enzymelist,kcatlist);
protein_efficiency_2 = yatp_2/totprotcost_2; %mmolATP/gProtein/h

% Arginine catabolism
[protcost_3,totprotcost_3] = CalculateProteinCost(rxnlist_3,fluxlist_3,Info_enzyme,kcat_glc,enzymelist,kcatlist);
protein_efficiency_3 = yatp_3/totprotcost_3; %mmolATP/gProtein/h


%% Small model simulation

% Estimate ATP consumption for growth, which is the sum of ATP for GAM,
% R_M_PROTS_LLA_v3, R_M_RNAS_LLA, and R_M_DNAS_LLA. The others can be
% neglected. 
GAM_MModel = 33.783;
newGAM = GAM_MModel + 18.0895 + 0.1316 + 0.10138;
NGAM_MModel = 3.72;

% Estimate proteome allocation to the three pathways.
tot_rxns = unique([rxnlist_1;rxnlist_2;rxnlist_3]);
tot_genes = cell(0,1);
for i = 1:length(tot_rxns)
    z = model.grRules{ismember(model.rxns,tot_rxns(i))};
    if ~isempty(z)
        z = strrep(z,'(',''); z = strrep(z,')','');
        z = strrep(z,'or',' '); z = strrep(z,'and',' ');
        z = strtrim(z);
        z = strsplit(z);%split multiple genes
    end
    tot_genes = [tot_genes;z'];
end
tot_genes = unique(tot_genes);
clear tot_rxns i z;

tot_proteome = 0.046;


% Estimate upper limit on arg uptake. The small model does not account for
% biomass formation, so arginine is only used to produce ATP.
argUB = 1.6568 - 0.1764; % total arg consumption minus arg composition in biomass.

% Simulations
mu_list = 0.01:0.01:0.8;

flux = zeros(3,length(mu_list));
rel_flux = zeros(2,length(mu_list));
unused_prot = zeros(1,length(mu_list));

for i = 1:length(mu_list)
    mu = mu_list(i);

    A = [-yatp_1 -yatp_2 -yatp_3;totprotcost_1 totprotcost_2 totprotcost_3];
    
    b = [-newGAM*mu-NGAM_MModel tot_proteome];

    Aeq = [];
    beq = [];

    lb = [0 0 0];
    ub = [100 100 argUB*mu];

    f = [1 1 0];

    x = linprog(f,A,b,Aeq,beq,lb,ub);
    
    flux(:,i) = x;
    rel_flux(:,i) = [x(1)/(x(1)+x(2));x(2)/(x(1)+x(2))];
    used_prot = totprotcost_1*x(1)+totprotcost_2*x(2)+totprotcost_3*x(3);
    unused_prot(i) = tot_proteome - round(used_prot,6);
end

clear i mu A b Aeq beq lb ub f x factor used_prot;

%% Plot
color_ma = [55,126,184]/255;
color_la = [228,26,28]/255;
color_arg = [152,78,163]/255;


figure('Name','1');
box on;
hold on;
plot(mu_list,rel_flux(1,:),'-','LineWidth',1.5,'Color',color_ma);
plot(mu_list,rel_flux(2,:),'-','LineWidth',1.5,'Color',color_la);
ylabel('Glucose fraction','FontSize',7,'FontName','Helvetica');
xlim([0 0.8]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
legend({'Mixed acid','Lactic acid'},'FontSize',7,'FontName','Helvetica','location','west');
set(gcf,'position',[0 600 160 100]);
set(gca,'position',[0.17 0.33 0.35 0.35]);

figure('Name','2');
box on;
plot(mu_list,flux(3,:),'-','LineWidth',1.5,'Color',color_arg);
ylim([0 1.2]);
xlim([0 0.8]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Arginine uptake','FontSize',7,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
% legend({'Arginine'},'FontSize',6,'FontName','Helvetica','location','nw');
set(gcf,'position',[0 400 160 100]);
set(gca,'position',[0.17 0.33 0.35 0.35]);

figure('Name','3');
box on;
plot(mu_list,unused_prot,'-','LineWidth',1.5,'Color','k');

xlim([0 0.8]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Inactive proteome','FontSize',7,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[0 200 160 100]);
set(gca,'position',[0.17 0.33 0.35 0.35]);














