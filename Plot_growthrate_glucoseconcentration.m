%% Max growth rates for a range of glucose concentrations

% Timing: ~ 2200 s

% With the saturation saturation factor are performed.

% Simulated results will be saved in the folder 'Results'.

% Figures can be obtained by running the codes starting from line 136.

tic;
load('pcLactis_Model.mat');
model = pcLactis_Model;

%% Optimization setting.

rxnID = 'R_dummy_assumed_Monomer';
osenseStr = 'Maximize';

%% Parameters.
GAM = 42;%ATP coefficient in the new biomass equation.
NGAM = 2.5; %(mmol/gCDW/h)
f_unmodeled = 0.42; %proportion of unmodeled protein in total protein (g/g)

model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
[model,f] = ChangeUnmodeledProtein(model,f_unmodeled);

kcat_glc = 180;%kcat value of glucose transporter
Km = 21;%Km of glucose transporter, unit: uM (PMID: 30630406)
f_transporter = 0.01;%fraction of glucose transporter in total proteome

%% Data import.
load('Info_enzyme.mat');
load('Info_mRNA.mat');
load('Info_protein.mat');
load('Info_ribosome.mat');
load('Info_tRNA.mat');
[~, ~, exchange_raw] = xlsread('Exchange_reaction_setting.xlsx','MaxMu');
[~, ~, aa_raw] = xlsread('Exchange_reaction_setting.xlsx','AA_factors');

Exchange_AAs = aa_raw(2:end,1);
LBfactor_AAs = cell2mat(aa_raw(2:end,2));
clear aa_raw;

%% Set reaction rates.
% Block uptake direction of all the exchange reactions.
idx = cellfun(@(x) contains(x,'R_M_EX_'),model.rxns,'UniformOutput',false);
M_exchange_reactions = model.rxns(cell2mat(idx));
model = changeRxnBounds(model,M_exchange_reactions,0,'l');
clear idx M_exchange_reactions;
% Set bounds for some exchange reactions.
Exchange_reactions = exchange_raw(2:end,1);
LB = cell2mat(exchange_raw(2:end,2));
UB = cell2mat(exchange_raw(2:end,3));
model = changeRxnBounds(model,Exchange_reactions,LB,'l');
model = changeRxnBounds(model,Exchange_reactions,UB,'u');
clear exchange_raw Exchange_reactions LB UB;
% Block some reactions in the M model.
model = changeRxnBounds(model,'R_M_biomass_LLA',0,'b');
model = changeRxnBounds(model,'R_M_biomass_LLA_atpm',0,'b');
model = changeRxnBounds(model,'R_M_PROTS_LLA',0,'b');
model = changeRxnBounds(model,'R_M_PROTS_LLA_v2',0,'b');
model = changeRxnBounds(model,'R_M_PROTS_LLA_v3',0,'b');
model = changeRxnBounds(model,'R_M_MGt2pp_rvs',0,'b');%block infinite h[e]

% Block other glucose transporters
model = changeRxnBounds(model,'R_M_GLCpts_2',0,'b');
model = changeRxnBounds(model,'R_M_GLCpermease_fwd',0,'b');

%% Main part.

glc_list = [5 10 20 40 60 80 100 200];%unit: /uM

% with saturation factor
% obtain the global saturation factor
load('Egsf1_result.mat');
x = global_saturation_factor_list(:,1);
y = global_saturation_factor_list(:,2);
x = x(y ~= 1);
y = y(y ~= 1);
sf_coeff = x\y;
clear x y;

res = zeros(length(glc_list),2);

for i = 1:length(glc_list)
    
    glc_conc = glc_list(i);
    factor_glc = glc_conc / (glc_conc + Km);
    
    mu_low = 0;
    mu_high = 0.8;
    
    while mu_high-mu_low > 0.01
        mu_mid = (mu_low+mu_high)/2;
        disp(['Glucose concentration = ' num2str(glc_conc) ' uM; mu = ' num2str(mu_mid)]);
        model = changeRxnBounds(model,'R_biomass_dilution',mu_mid,'b');
        model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*mu_mid,'l');
        
        factor_k = sf_coeff * mu_mid;
        if factor_k > 1
            factor_k = 1;
        end
        
        fileName = WriteLPSatFactor(model,mu_mid,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,factor_glc,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);

        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,solME_status,~] = ReadSoplexResult(fileName_out,model);
        
        if strcmp(solME_status,'optimal')
            mu_low = mu_mid;
        else
            mu_high = mu_mid;
        end
    end
    
    res(i,:) = [glc_conc mu_low];
end

cd Results/;
save('Pgg_result.mat','res');
cd ../;

clear;
toc;

%% Figure.
load('Pgg_result.mat');

figure('Name','1');
hold on;
box on;
% plot simulations
x = res(:,1);
y = res(:,2);
plot(x,y,'o','MarkerSize',5,'Color',[178,24,43]/255);
% plot measurement-based equation (PMID: 30630406)
s = 0:1:200;
gr = 0.71*s./(s+19);
plot(s,gr,'-','LineWidth',1,'Color',[77,77,77]/255);
set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Glucose concentration (uM)','FontSize',14,'FontName','Helvetica');
ylabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');

ylim([0 0.8]);

legend({'Simulated','Experimental'},'FontSize',12,'FontName','Helvetica','location','se');

set(gcf,'position',[200 0 230 130]);
set(gca,'position',[0.23 0.3 0.72 0.64]);

clear x y s gr;

