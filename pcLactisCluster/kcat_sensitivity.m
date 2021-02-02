%% Sensitivity analysis for kcats of genome-scale model.
function kcat_sensitivity(fold,glc_conc,a,b)

addpath(genpath('../pcLactisCluster/')); %add to path
addpath(genpath('/proj/nobackup/snic2021-22-21/tools/'));
workdir = pwd;
cd '/proj/nobackup/snic2021-22-21/tools/cobratoolbox';
initCobraToolbox;
savepath '~/pathdef.m';
cd(workdir);

soplexlocstr = '/proj/nobackup/snic2021-22-21/tools/soplex-4.0.0/build/bin/soplex';
soplexsettingstr = ' -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s';

load('pcLactis_Model_Cluster.mat');
model = pcLactis_Model;

load('param_list.mat');

%% Optimization setting.

rxnID = 'R_dummy_assumed_Monomer';
osenseStr = 'Maximize';

%% Parameters.
GAM = 36;%ATP coefficient in the new biomass equation.
NGAM = 2; %(mmol/gCDW/h)

f_unmodeled = 0.4;
f_transporter = 0.0082;
[model,f] = ChangeUnmodeledProteinCluster(model,f_unmodeled);

model = ChangeATPinBiomassCluster(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');

kcat_glc = 180;%kcat value of glucose transporter
Km = 21;%Km of glucose transporter, unit: uM (PMID: 30630406)

%% Data import.
load('Info_enzyme_Cluster.mat');
load('Info_mRNA_Cluster.mat');
load('Info_protein_Cluster.mat');
load('Info_ribosome_Cluster.mat');
[~, ~, exchange_raw] = xlsread('Exchange_reaction_setting_Cluster.xlsx','MaxMu');
[~, ~, aa_raw] = xlsread('Exchange_reaction_setting_Cluster.xlsx','AA_factors');

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
model = changeRxnBounds(model,'R_M_PROTS_LLA_v3',0,'b');

% Block other glucose transporters
model = changeRxnBounds(model,'R_M_GLCpts_1',0,'b');
model = changeRxnBounds(model,'R_M_GLCt2_fwd',0,'b');

% Block one of ADH isozymes llmg_0955
model = changeRxnBounds(model,'R_M_ALCD2x_1_rvs',0,'b');

% Block pyruvate oxidase
model = changeRxnBounds(model,'R_M_PYROX_1',0,'b');

%% Main part.
precision = 1e-6;
org_min_mu = 0;
org_max_mu = 1.1;
factor_k = 1;

factor_glc = glc_conc / (glc_conc + Km);

if fold > 1
    direction = 'up';
elseif fold < 1
    direction = 'dn';
end

if b >= length(param_list)
    b = length(param_list);
end

mu_list = zeros(b-a+1,1);
for i = a:b
    ID = param_list{i};
    disp([num2str(i) '/' num2str(length(param_list)) ': ' ID]);
    model_tmp = model;
    
    label_tmp = [direction '_' num2str(glc_conc) '_' num2str(i)];
    
    mu_low = org_min_mu;
    mu_high = org_max_mu;
    
    while mu_high-mu_low > precision
        mu_mid = (mu_low+mu_high)/2;
        disp(['mu = ' num2str(mu_mid)]);
        model_tmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_mid,'b');
        model_tmp = changeRxnBounds(model_tmp,Exchange_AAs,LBfactor_AAs*mu_mid,'l');
        fileName = WriteLPkcatSensitivityCluster(model_tmp,mu_mid,f,osenseStr,rxnID,factor_k,...
                                                 f_transporter,kcat_glc,factor_glc,...
                                                 Info_enzyme,...
                                                 Info_mRNA,...
                                                 Info_protein,...
                                                 Info_ribosome,...
                                                 ID,fold,label_tmp);
        command = sprintf([soplexlocstr soplexsettingstr],fileName,fileName);
        system(command,'-echo');
        fileName_out = strcat('Simulation_',label_tmp,'.lp.out');
        [~,solME_status,solME_full] = ReadSoplexResultCluster(fileName_out,model_tmp);
        if strcmp(solME_status,'optimal')
            mu_low = mu_mid;
            flux_tmp = solME_full;
        else
            mu_high = mu_mid;
        end
    end
    mu_list(i-a+1,1) = mu_low;
end

filename = strcat('kcat_sens_',direction,'_',num2str(glc_conc),'_',num2str(a),'_',num2str(b),'.mat');
cd ClusterResults/;
save(filename,'mu_list');
cd ../;

end
