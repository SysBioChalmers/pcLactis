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
clear num txt;

% Glycolysis + mixed acid fermentation
rxns_tmp = rxnlist_1;
flux_tmp = fluxlist_1;
protcost_1 = zeros(length(rxns_tmp),1);
for i = 1:length(rxns_tmp)
    rxnid = rxns_tmp{i};
    ezmid = strcat(rxnid(3:end),'_');
    flux = flux_tmp(i);
    if ismember({rxnid},{'R_M_GLCpts'})
        kcat_tmp = kcat_glc * 3600;%/h
        mw_tmp = Info_enzyme.MW(ismember(Info_enzyme.ID,'M_GLCpts_2_Enzyme_c'));
        protcost_tmp = mw_tmp/kcat_tmp;
    else
        if any(contains(enzymelist,ezmid))
            enzymes_tmp = enzymelist(contains(enzymelist,ezmid));
            kcats_tmp = kcatlist(contains(enzymelist,ezmid));
            if flux < 0
                idx_tmp = contains(enzymes_tmp,'_rvs');
            else
                idx_tmp = ~contains(enzymes_tmp,'_rvs');
            end
            enzymes_tmptmp = enzymes_tmp(idx_tmp);
            kcats_tmptmp = kcats_tmp(idx_tmp);
            [~, b] = ismember(enzymes_tmptmp,Info_enzyme.ID);
            mws_tmp = Info_enzyme.MW(b);
            pcs_tmp = mws_tmp./kcats_tmptmp;
            protcost_tmp = min(pcs_tmp);
        else
            protcost_tmp = 0;
        end
    end
    protcost_1(i) = protcost_tmp; %gprotein/gCDW per flux(mol/gCDW/h)
end
totprotcost_1 = sum(protcost_1.*abs(fluxlist_1))/1000; %gprotein/gCDW per flux(mmol/gCDW/h) of glucose
protein_efficiency_1 = yatp_1/totprotcost_1; %mmolATP/gProtein/h

% Glycolysis + lactic acid fermentation
rxns_tmp = rxnlist_2;
flux_tmp = fluxlist_2;
protcost_2 = zeros(length(rxns_tmp),1);
for i = 1:length(rxns_tmp)
    rxnid = rxns_tmp{i};
    ezmid = strcat(rxnid(3:end),'_');
    flux = flux_tmp(i);
    if ismember({rxnid},{'R_M_GLCpts'})
        kcat_tmp = kcat_glc * 3600;%/h
        mw_tmp = Info_enzyme.MW(ismember(Info_enzyme.ID,'M_GLCpts_2_Enzyme_c'));
        protcost_tmp = mw_tmp/kcat_tmp;
    else
        if any(contains(enzymelist,ezmid))
            enzymes_tmp = enzymelist(contains(enzymelist,ezmid));
            kcats_tmp = kcatlist(contains(enzymelist,ezmid));
            if flux < 0
                idx_tmp = contains(enzymes_tmp,'_rvs');
            else
                idx_tmp = ~contains(enzymes_tmp,'_rvs');
            end
            enzymes_tmptmp = enzymes_tmp(idx_tmp);
            kcats_tmptmp = kcats_tmp(idx_tmp);
            [~, b] = ismember(enzymes_tmptmp,Info_enzyme.ID);
            mws_tmp = Info_enzyme.MW(b);
            pcs_tmp = mws_tmp./kcats_tmptmp;
            protcost_tmp = min(pcs_tmp);
        else
            protcost_tmp = 0;
        end
    end
    protcost_2(i) = protcost_tmp; %gprotein/gCDW per flux(mol/gCDW/h)
end
totprotcost_2 = sum(protcost_2.*abs(fluxlist_2))/1000; %gprotein/gCDW per flux(mmol/gCDW/h) of glucose
protein_efficiency_2 = yatp_2/totprotcost_2; %mmolATP/gProtein/h

% Arginine catabolism
rxns_tmp = rxnlist_3;
flux_tmp = fluxlist_3;
protcost_3 = zeros(length(rxns_tmp),1);
for i = 1:length(rxns_tmp)
    rxnid = rxns_tmp{i};
    ezmid = strcat(rxnid(3:end),'_');
    flux = flux_tmp(i);
    if ismember({rxnid},{'R_M_GLCpts'})
        kcat_tmp = kcat_glc * 3600;%/h
        mw_tmp = Info_enzyme.MW(ismember(Info_enzyme.ID,'M_GLCpts_2_Enzyme_c'));
        protcost_tmp = mw_tmp/kcat_tmp;
    else
        if any(contains(enzymelist,ezmid))
            enzymes_tmp = enzymelist(contains(enzymelist,ezmid));
            kcats_tmp = kcatlist(contains(enzymelist,ezmid));
            if flux < 0
                idx_tmp = contains(enzymes_tmp,'_rvs');
            else
                idx_tmp = ~contains(enzymes_tmp,'_rvs');
            end
            enzymes_tmptmp = enzymes_tmp(idx_tmp);
            kcats_tmptmp = kcats_tmp(idx_tmp);
            [~, b] = ismember(enzymes_tmptmp,Info_enzyme.ID);
            mws_tmp = Info_enzyme.MW(b);
            pcs_tmp = mws_tmp./kcats_tmptmp;
            protcost_tmp = min(pcs_tmp);
        else
            protcost_tmp = 0;
        end
    end
    protcost_3(i) = protcost_tmp; %gprotein/gCDW per flux(mol/gCDW/h)
end
totprotcost_3 = sum(protcost_3.*abs(fluxlist_3))/1000; %gprotein/gCDW per flux(mmol/gCDW/h) of arginine
protein_efficiency_3 = yatp_3/totprotcost_3; %mmolATP/gProtein/h

clear b enzymelist enzymelist_tmp enzymes_tmp enzymes_tmptmp;
clear ezmid flux flux_tmp i idx_tmp kcat_glc kcat_tmp kcatlist kcatlist_tmp;
clear kcats_tmp kcats_tmptmp mw_tmp mws_tmp pcs_tmp protcost_tmp rxnid;
clear rxns_tmp;

%% Small model simulation
mu_list = 0.1:0.025:0.725;

flux = [];
rel_flux = [];

for i = 1:length(mu_list)
    mu = mu_list(i);

    A = [-yatp_1 -yatp_2 -yatp_3;totprotcost_1 totprotcost_2 totprotcost_3];
    
    b = [-80.44*mu+1.64 0.06];

    Aeq = [];
    beq = [];

    lb = [0 0 0];
    ub = [100 100 1.6568*mu];

    f = [1 1 0];

    x = linprog(f,A,b,Aeq,beq,lb,ub);
    
    flux = [flux x];
    rel_flux = [rel_flux [x(1)/(x(1)+x(2));x(2)/(x(1)+x(2))]];
    
end

clear i mu A b Aeq beq lb ub f x factor;

figure('Name','1');
hold on;
plot(mu_list,rel_flux(1,:),'-o','LineWidth',1.5,'Color',[55,126,184]/256,'MarkerSize',8,'MarkerEdgeColor',[55,126,184]/256,'MarkerFaceColor',[55,126,184]/256);
plot(mu_list,rel_flux(2,:),'-o','LineWidth',1.5,'Color',[228,26,28]/256,'MarkerSize',8,'MarkerEdgeColor',[228,26,28]/256,'MarkerFaceColor',[228,26,28]/256);
ylabel('Fraction in total glucose carbon','FontSize',12,'FontName','Helvetica');

set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',12,'FontName','Helvetica');
legend({'Mixed acid','Lactate'},'FontSize',12,'FontName','Helvetica','location','west');
set(gcf,'position',[0 0 335 300]);
set(gca,'position',[0.11 0.11 0.75 0.8]);

figure('Name','2');
plot(mu_list,flux(3,:),'-o','LineWidth',1.5,'Color',[77,175,74]/256,'MarkerSize',8,'MarkerEdgeColor',[77,175,74]/256,'MarkerFaceColor',[77,175,74]/256);
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('q(mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',12,'FontName','Helvetica');
legend({'Arginine'},'FontSize',12,'FontName','Helvetica','location','nw');
set(gcf,'position',[0 0 335 300]);
set(gca,'position',[0.12 0.11 0.75 0.8]);
















