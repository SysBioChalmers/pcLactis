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
rxnlist_1 = model_1.rxns(abs(sol_1.x) > 1e-5);
fluxlist_1 = sol_1.x(abs(sol_1.x) > 1e-5);
yatp_1 = sol_1.x(ismember(model_1.rxns,'R_M_ATPM')) / -sol_1.x(ismember(model_1.rxns,'R_M_EX_glc__D_e'));
clear model_1 sol_1;

%% Glycolysis + lactic acid fermentation
model_2 = changeObjective(model,'R_M_EX_lac__L_e',1);
model_2 = changeRxnBounds(model_2,'R_M_EX_glc__D_e',-1,'b');
model_2 = changeRxnBounds(model_2,'R_M_ATPM',2,'b');
sol_2 = optimizeCbModel(model_2,'max','one');
rxnlist_2 = model_2.rxns(abs(sol_2.x) > 1e-5);
fluxlist_2 = sol_2.x(abs(sol_2.x) > 1e-5);
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
rxnlist_3 = model_3.rxns(abs(sol_3.x) > 1e-5);
fluxlist_3 = sol_3.x(abs(sol_3.x) > 1e-5);
yatp_3 = sol_3.x(ismember(model_3.rxns,'R_M_ATPM')) / -sol_3.x(ismember(model_3.rxns,'R_M_EX_arg__L_e'));
% clear model_3 sol_3;

%% Minimal protein cost for each rxn and pathway
load('Info_enzyme.mat');

load('kcat_M.mat');
enzymelist_tmp = kcat_M.gpr;
kcatlist_tmp = kcat_M.kcat;

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

% value = yatp_3-totprotcost_3/(totprotcost_1-totprotcost_2);

%% Small model simulation

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

tot_proteome = 0.46*0.1;

% Fix glucose uptake maximize ATP production
% Estimate upper limit on arg uptake. The small model does not account for
% biomass formation, so arginine is only used to produce ATP.
argUB = 0.0512; % total arg consumption minus arg composition in biomass.

qglc_list = 0.1:0.1:20.6;

%% Simulations Ref
flux = zeros(4,length(qglc_list));
unused_prot = zeros(1,length(qglc_list));
f = [-yatp_1 -yatp_2 -yatp_3];
A = [1 1 0;totprotcost_1 totprotcost_2 totprotcost_3];
Aeq = [];
beq = [];
lb = [0 0 0];
for i = 1:length(qglc_list)
    qglc = qglc_list(i);
    
    b = [qglc tot_proteome];
    
    ub = [100 100 argUB*qglc];

    [x, fval] = linprog(f,A,b,Aeq,beq,lb,ub);
    
    flux(:,i) = [x;-fval];
    used_prot = totprotcost_1*x(1)+totprotcost_2*x(2)+totprotcost_3*x(3);
    unused_prot(i) = tot_proteome - round(used_prot,6);
end

increasefactor = 0.01;

% 1. increase bound of arg uptake
res_arg = zeros(5,length(qglc_list)); %[ref_arg;new_arg;ref_qatp;new_qatp;qs]

for i = 1:length(qglc_list)
    b = [qglc_list(i) tot_proteome];
    ub = [100 100 argUB*qglc_list(i)*(1+increasefactor)];
    [x, fval] = linprog(f,A,b,Aeq,beq,lb,ub);
    res_arg(1,i) = argUB*qglc_list(i);
    res_arg(2,i) = argUB*qglc_list(i)*(1+increasefactor);
    res_arg(3,i) = flux(4,i);
    res_arg(4,i) = -fval;
    res_arg(5,i) = qglc_list(i);
end
sensitivity_arg = (res_arg(4,:)-res_arg(3,:))./(res_arg(2,:)-res_arg(1,:));
sensitivity_arg = round(sensitivity_arg,6);
sensitivity_arg = sensitivity_arg.*round((res_arg(1,:)./res_arg(3,:)),6);

% 2. increase bound of glucose uptake, similar to increase glucose
% transporter.
res_glc = zeros(5,length(qglc_list)); %[ref_glc;new_glc;ref_qatp;new_qatp;qs]

for i = 1:length(qglc_list)
    b = [qglc_list(i)*(1+increasefactor) tot_proteome];
    ub = [100 100 argUB*qglc_list(i)];
    [x, fval] = linprog(f,A,b,Aeq,beq,lb,ub);
    res_glc(1,i) = qglc_list(i);
    res_glc(2,i) = qglc_list(i)*(1+increasefactor);
    res_glc(3,i) = flux(4,i);
    res_glc(4,i) = -fval;
    res_glc(5,i) = qglc_list(i);
end
sensitivity_glc = (res_glc(4,:)-res_glc(3,:))./(res_glc(2,:)-res_glc(1,:));
sensitivity_glc = round(sensitivity_glc,6);
sensitivity_glc = sensitivity_glc.*round((res_glc(1,:)./res_glc(3,:)),6);

% 3. increase constraint on proteome
res_proteome = zeros(5,length(qglc_list)); %[ref_proteome;new_proteome;ref_qatp;new_qatp;qs]

for i = 1:length(qglc_list)
    b = [qglc_list(i) tot_proteome*(1+increasefactor)];
    ub = [100 100 argUB*qglc_list(i)];
    [x, fval] = linprog(f,A,b,Aeq,beq,lb,ub);
    res_proteome(1,i) = tot_proteome;
    res_proteome(2,i) = totprotcost_1*x(1)+totprotcost_2*x(2)+totprotcost_3*x(3);
    res_proteome(3,i) = flux(4,i);
    res_proteome(4,i) = -fval;
    res_proteome(5,i) = qglc_list(i);
end
sensitivity_proteome = (res_proteome(4,:)-res_proteome(3,:))./(res_proteome(2,:)-res_proteome(1,:));
sensitivity_proteome = round(sensitivity_proteome,6);
sensitivity_proteome = sensitivity_proteome.*round((res_proteome(1,:)./res_proteome(3,:)),6);

%% Simulations new
% sensitivity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_totprotcost_1 = totprotcost_1 * 1;
new_totprotcost_2 = totprotcost_2 * 1;
new_totprotcost_3 = totprotcost_3 * 1;
new_tot_proteome = tot_proteome * 1;
new_argUB = argUB * 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

new_flux = zeros(4,length(qglc_list));
new_unused_prot = zeros(1,length(qglc_list));
f = [-yatp_1 -yatp_2 -yatp_3];
A = [1 1 0;new_totprotcost_1 new_totprotcost_2 new_totprotcost_3];
Aeq = [];
beq = [];
lb = [0 0 0];
for i = 1:length(qglc_list)
    qglc = qglc_list(i);
    
    b = [qglc new_tot_proteome];
    
    ub = [100 100 new_argUB*qglc];

    [x, fval] = linprog(f,A,b,Aeq,beq,lb,ub);
    
    new_flux(:,i) = [x;-fval];
    used_prot = new_totprotcost_1*x(1)+new_totprotcost_2*x(2)+new_totprotcost_3*x(3);
    new_unused_prot(i) = new_tot_proteome - round(used_prot,6);
end

increasefactor = 0.01;

% 1. increase bound of arg uptake
new_res_arg = zeros(5,length(qglc_list)); %[ref_arg;new_arg;ref_qatp;new_qatp;qs]

for i = 1:length(qglc_list)
    b = [qglc_list(i) new_tot_proteome];
    ub = [100 100 new_argUB*qglc_list(i)*(1+increasefactor)];
    [x, fval] = linprog(f,A,b,Aeq,beq,lb,ub);
    new_res_arg(1,i) = new_argUB*qglc_list(i);
    new_res_arg(2,i) = new_argUB*qglc_list(i)*(1+increasefactor);
    new_res_arg(3,i) = new_flux(4,i);
    new_res_arg(4,i) = -fval;
    new_res_arg(5,i) = qglc_list(i);
end
new_sensitivity_arg = (new_res_arg(4,:)-new_res_arg(3,:))./(new_res_arg(2,:)-new_res_arg(1,:));
new_sensitivity_arg = round(new_sensitivity_arg,6);
new_sensitivity_arg = new_sensitivity_arg.*round((new_res_arg(1,:)./new_res_arg(3,:)),6);

% 2. increase bound of glucose uptake, similar to increase glucose
% transporter.
new_res_glc = zeros(5,length(qglc_list)); %[ref_glc;new_glc;ref_qatp;new_qatp;qs]

for i = 1:length(qglc_list)
    b = [qglc_list(i)*(1+increasefactor) new_tot_proteome];
    ub = [100 100 new_argUB*qglc_list(i)];
    [x, fval] = linprog(f,A,b,Aeq,beq,lb,ub);
    new_res_glc(1,i) = qglc_list(i);
    new_res_glc(2,i) = qglc_list(i)*(1+increasefactor);
    new_res_glc(3,i) = new_flux(4,i);
    new_res_glc(4,i) = -fval;
    new_res_glc(5,i) = qglc_list(i);
end
new_sensitivity_glc = (new_res_glc(4,:)-new_res_glc(3,:))./(new_res_glc(2,:)-new_res_glc(1,:));
new_sensitivity_glc = round(new_sensitivity_glc,6);
new_sensitivity_glc = new_sensitivity_glc.*round((new_res_glc(1,:)./new_res_glc(3,:)),6);

% 3. increase constraint on proteome
new_res_proteome = zeros(5,length(qglc_list)); %[ref_proteome;new_proteome;ref_qatp;new_qatp;qs]

for i = 1:length(qglc_list)
    b = [qglc_list(i) new_tot_proteome*(1+increasefactor)];
    ub = [100 100 new_argUB*qglc_list(i)];
    [x, fval] = linprog(f,A,b,Aeq,beq,lb,ub);
    new_res_proteome(1,i) = new_tot_proteome;
    new_res_proteome(2,i) = new_totprotcost_1*x(1)+new_totprotcost_2*x(2)+new_totprotcost_3*x(3);
    new_res_proteome(3,i) = new_flux(4,i);
    new_res_proteome(4,i) = -fval;
    new_res_proteome(5,i) = qglc_list(i);
end
new_sensitivity_proteome = (new_res_proteome(4,:)-new_res_proteome(3,:))./(new_res_proteome(2,:)-new_res_proteome(1,:));
new_sensitivity_proteome = round(new_sensitivity_proteome,6);
new_sensitivity_proteome = new_sensitivity_proteome.*round((new_res_proteome(1,:)./new_res_proteome(3,:)),6);

%% Plot
color_ma = [55,126,184]/255;
color_la = [228,26,28]/255;
color_arg = [152,78,163]/255;


figure('Name','1');

subplot(6,1,1);
box on;
hold on;
plot(qglc_list,new_flux(1,:),'-','LineWidth',1,'Color',color_ma);
plot(qglc_list,new_flux(2,:),'-','LineWidth',1,'Color',color_la);
plot(qglc_list,flux(1,:),'-','LineWidth',1,'Color',[color_ma,0.2]);
plot(qglc_list,flux(2,:),'-','LineWidth',1,'Color',[color_la,0.2]);
ylim([0 25]);
xlim([0 21]);
set(gca, 'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');

subplot(6,1,2);
box on;
hold on;
plot(qglc_list,new_flux(3,:),'-','LineWidth',1,'Color',color_arg);
plot(qglc_list,flux(3,:),'-','LineWidth',1,'Color',[color_arg,0.2]);
ylabel('Flux','FontSize',6,...
       'FontName','Helvetica','Color','k');
set(gca, 'YColor','k');
ylim([0 0.99]);
xlim([0 21]);
set(gca,'FontSize',6,'FontName','Helvetica');

subplot(6,1,3);
box on;
hold on;
plot(qglc_list,new_unused_prot,'-','LineWidth',1,'Color',[0,0,0]);
plot(qglc_list,unused_prot,'-','LineWidth',1,'Color',[0,0,0,0.2]);
ylabel('Content','FontSize',6,...
       'FontName','Helvetica','Color','k');
xlim([0 21]);
ylim([0 0.099]);
set(gca,'FontSize',6,'FontName','Helvetica');
set(gca, 'YColor','k');
set(gcf,'position',[400 500 80 150]);

subplot(6,1,4);
box on;
hold on;
plot(new_res_glc(5,:),new_sensitivity_glc,'-','LineWidth',1,'Color',[197,27,138]/255);
plot(new_res_glc(5,:),sensitivity_glc,'-','LineWidth',1,'Color',[[197,27,138]/255,0.2]);
xlim([0 21]);
ylim([0 1.2]);
set(gca, 'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
title('Glucose uptake','FontSize',6,'FontName','Helvetica');

subplot(6,1,5);
box on;
hold on;
plot(new_res_proteome(5,:),new_sensitivity_proteome,'-','LineWidth',1,'Color',[197,27,138]/255);
plot(new_res_proteome(5,:),sensitivity_proteome,'-','LineWidth',1,'Color',[[197,27,138]/255,0.2]);
xlim([0 21]);
ylim([0 1.2]);
ylabel('Sensitivity score','FontSize',6,...
       'FontName','Helvetica','Color','k');
set(gca,'FontSize',6,'FontName','Helvetica');
title('Total proteome','FontSize',6,'FontName','Helvetica');
set(gca, 'YColor','k');

subplot(6,1,6);
box on;
hold on;
plot(new_res_arg(5,:),new_sensitivity_arg,'-','LineWidth',1,'Color',[197,27,138]/255);
plot(new_res_arg(5,:),sensitivity_arg,'-','LineWidth',1,'Color',[[197,27,138]/255,0.2]);
xlim([0 21]);
ylim([0 0.06]);
set(gca, 'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
title('Arginine uptake','FontSize',6,'FontName','Helvetica');
xlabel('Glucose uptake','FontSize',6,'FontName','Helvetica','Color','k');
set(gca, 'XColor','k');
set(gcf,'position',[400 200 80 300]);

