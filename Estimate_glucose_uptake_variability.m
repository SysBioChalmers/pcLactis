%% FVA for glucose uptake rate 
% Calculate UB and LB of the glucose uptake rate for each glucose
% concentration.

% Timing: ~ 400 s

% Simulated results will be saved in the folder 'Results'.

% Figures can be obtained by running the codes starting from line 149.

tic;
load('pcLactis_Model.mat');
model = pcLactis_Model;

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

% obtain the global saturation factor
load('Egsf2_result.mat');
x = global_saturation_factor_list(:,1);
y = global_saturation_factor_list(:,2);
x = x(y ~= 1);
y = y(y ~= 1);
sf_coeff = x\y;
clear x y;

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

%% Main part.
load('Sglc_result_with_sf.mat');

D_list = glc_conc_with_sf(:,1);%unit: /h
glc_conc_list = glc_conc_with_sf(:,2);%unit: uM

fluxes_simulated = zeros(length(model.rxns),2*length(D_list));

for i = 1:length(D_list)
    
    D = D_list(i);
	factor_k = sf_coeff * D;
    if factor_k > 1
        factor_k = 1;
    end
	model = changeRxnBounds(model,'R_biomass_dilution',D,'b');
	model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*D,'l');
    
	glc_conc = glc_conc_list(i);
    factor_glc = glc_conc / (glc_conc + Km);
    
    rxnID = 'R_M_EX_glc_LPAREN_e_RPAREN_'; 
    
    % Minimize glucose uptake rate
    osenseStr = 'Maximize';
    disp(['D = ' num2str(D) '; Minimizing glucose uptake']);
    
	fileName = WriteLPSatFactor(model,D,f,osenseStr,rxnID,factor_k,...
                                f_transporter,kcat_glc,factor_glc,...
                                Info_enzyme,...
                                Info_mRNA,...
                                Info_protein,...
                                Info_ribosome,...
                                Info_tRNA);
                   
	command =sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
	system(command,'-echo');
	fileName_out = 'Simulation.lp.out';
	[~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model);
    if strcmp(solME_status,'optimal')
        fluxes_simulated(:,2*i-1) = solME_full;
    else
        fluxes_simulated(:,2*i-1) = zeros(length(model.rxns),1);
    end
    
    % Maximize glucose uptake rate 
    osenseStr = 'Minimize';
    disp(['D = ' num2str(D) '; Maximizing glucose uptake']);
    
	fileName = WriteLPSatFactor(model,D,f,osenseStr,rxnID,factor_k,...
                                f_transporter,kcat_glc,factor_glc,...
                                Info_enzyme,...
                                Info_mRNA,...
                                Info_protein,...
                                Info_ribosome,...
                                Info_tRNA);
                   
	command =sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
	system(command,'-echo');
	fileName_out = 'Simulation.lp.out';
	[~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model);
    
    if strcmp(solME_status,'optimal')
        fluxes_simulated(:,2*i) = solME_full;
    else
        fluxes_simulated(:,2*i) = zeros(length(model.rxns),1);
    end
end

cd Results/;
save('Eguv_fluxes.mat','fluxes_simulated');
cd ../;

clear;

toc;

%% Figure.
load('Eguv_fluxes.mat');
load('pcLactis_Model.mat');
model = pcLactis_Model;

mu = fluxes_simulated(strcmp(model.rxns,'R_biomass_dilution'),:);
glc = -fluxes_simulated(strcmp(model.rxns,'R_M_EX_glc_LPAREN_e_RPAREN_'),:);
mu_list = mu(1,1:2:end);
low_glc_list = glc(1,1:2:end);
high_glc_list = glc(1,2:2:end);
aver_glc_list = (low_glc_list + high_glc_list)/2;
rel_aver_list = 100 * ones(length(aver_glc_list),1);
rel_err_list = 100 * (high_glc_list - low_glc_list)/2 ./ aver_glc_list;

% Relative values
figure('Name','1');
errorbar(mu_list,rel_aver_list,rel_err_list,'.','LineWidth',1);

set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel('Relative q_s (%)','FontSize',14,'FontName','Helvetica');

ylim([98.5,101.5]);

set(gcf,'position',[0 300 230 130]);
set(gca,'position',[0.27 0.29 0.68 0.66]);

% Absolute values
figure('Name','2');
hold on;
box on;
plot(mu_list,aver_glc_list,'o','LineWidth',1);
set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel(['Absolute q_s',char(13,10)','(mmol/gCDW/h)'],'FontSize',14,'FontName','Helvetica');

set(gcf,'position',[0 0 230 130]);
set(gca,'position',[0.27 0.29 0.68 0.66]);

clear;
