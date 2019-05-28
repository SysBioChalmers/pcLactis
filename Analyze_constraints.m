%% Sensitivity analysis for total modeled proteome and glucose transporter.

% Timing: ~ 60000 s

% With the saturation saturation factor are performed.

% Simulated results will be saved in the folder 'Results'.

% Figures can be obtained by running the codes starting from line .

tic;
load('pcLactis_Model.mat');
model = pcLactis_Model;

%% Optimization setting.

rxnID = 'R_dummy_assumed_Monomer';
osenseStr = 'Maximize';

%% Parameters.
GAM = 42;%ATP coefficient in the new biomass equation.
NGAM = 2.5; %(mmol/gCDW/h)

model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');

kcat_glc = 180;%kcat value of glucose transporter
Km = 21;%Km of glucose transporter, unit: uM (PMID: 30630406)

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

glc_list = [2 4 6 8 10 20 40 60 80 100 200 1000 10000 100000 1000000];%unit: /uM
modeled_protein_list = (1-0.42)*(1:0.1:1.5);
glucose_transporter_list = 0.01*(1:0.1:1.5);
res1 = zeros(length(modeled_protein_list),4,length(glc_list));
res2 = zeros(length(glucose_transporter_list),4,length(glc_list));

% obtain the global saturation factor
load('Egsf1_result.mat');
x = global_saturation_factor_list(:,1);
y = global_saturation_factor_list(:,2);
x = x(y ~= 1);
y = y(y ~= 1);
sf_coeff = x\y;
clear x y;

for i = 1:length(glc_list)
    glc_conc = glc_list(i);
    factor_glc = glc_conc / (glc_conc + Km);
    
    for j = 1:length(modeled_protein_list)
        f_unmodeled = 1-modeled_protein_list(j);
        [model,f] = ChangeUnmodeledProtein(model,f_unmodeled);
        f_transporter = 0.01;
        
        mu_low = 0;
        mu_high = 2;

        while mu_high-mu_low > 0.01
            mu_mid = (mu_low+mu_high)/2;
            disp(['Glucose concentration = ' num2str(glc_conc) ' uM; modeled protein = ' num2str((1-f_unmodeled)) '; glucose transporter = ' num2str(f_transporter) '; mu = ' num2str(mu_mid)]);
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
        
        model = changeRxnBounds(model,'R_biomass_dilution',mu_low,'b');
        model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*mu_low,'l');
        factor_k = sf_coeff * mu_low;
        if factor_k > 1
            factor_k = 1;
        end
        fileName = WriteLPSatFactor(model,mu_low,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,factor_glc,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,~,solME_full] = ReadSoplexResult(fileName_out,model);
        
        mu = solME_full(strcmp(model.rxns,'R_biomass_dilution'),1);
        glc = -solME_full(strcmp(model.rxns,'R_M_EX_glc_LPAREN_e_RPAREN_'),1);
        arg = -solME_full(strcmp(model.rxns,'R_M_EX_arg_L_LPAREN_e_RPAREN_'),1);
        res1(j,:,i) = [1-f_unmodeled mu mu/(glc*180/1000) arg*174/1000/mu];
                                       %g_CDW/g_glucose   g_arginine/g_CDW
    end
        
	for k = 1:length(glucose_transporter_list)
        f_transporter = glucose_transporter_list(k);
        f_unmodeled = 0.42;
        [model,f] = ChangeUnmodeledProtein(model,f_unmodeled);
        
        mu_low = 0;
        mu_high = 2;
        
        while mu_high-mu_low > 0.01
            mu_mid = (mu_low+mu_high)/2;
            disp(['Glucose concentration = ' num2str(glc_conc) ' uM; modeled protein = ' num2str((1-f_unmodeled)) '; glucose transporter = ' num2str(f_transporter) '; mu = ' num2str(mu_mid)]);
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
        
        model = changeRxnBounds(model,'R_biomass_dilution',mu_low,'b');
        model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*mu_low,'l');
        factor_k = sf_coeff * mu_low;
        if factor_k > 1
            factor_k = 1;
        end
        fileName = WriteLPSatFactor(model,mu_low,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,factor_glc,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,~,solME_full] = ReadSoplexResult(fileName_out,model);
        
        mu = solME_full(strcmp(model.rxns,'R_biomass_dilution'),1);
        glc = -solME_full(strcmp(model.rxns,'R_M_EX_glc_LPAREN_e_RPAREN_'),1);
        arg = -solME_full(strcmp(model.rxns,'R_M_EX_arg_L_LPAREN_e_RPAREN_'),1);
        res2(k,:,i) = [f_transporter mu mu/(glc*180/1000) arg*174/1000/mu];
                                       %g_CDW/g_glucose   g_arginine/g_CDW
	end
end
    


cd Results/;
save('Ac_result1.mat','res1');
save('Ac_result2.mat','res2');
cd ../;

clear;
toc;

%% Figures
load('Ac_result1.mat');
load('Ac_result2.mat');
glc_list = [2 4 6 8 10 20 40 60 80 100 200 1000 10000 100000 1000000];%unit: /uM
idx_1 = 1:5;
idx_2 = 6:11;
idx_3 = 12:length(glc_list);
st1 = glc_list(idx_1);%unit: /uM
st2 = glc_list(idx_2);%unit: /uM
st3 = glc_list(idx_3);%unit: /uM

clr_mu = [31,120,180]/255;
clr_BiomassGlc = [51,160,44]/255;
clr_ArgBiomass = [227,26,28]/255;

% Reference
ref_mu = permute(res1(1,2,:),[3 2 1]);
ref_BiomassGlc = permute(res1(1,3,:),[3 2 1]);
ref_ArgBiomass = permute(res1(1,4,:),[3 2 1]);

figure();
hold on;
plot(glc_list(idx_1),ref_mu(idx_1),'Color',clr_mu);
plot(glc_list(idx_1),ref_BiomassGlc(idx_1),'Color',clr_BiomassGlc);
plot(glc_list(idx_1),ref_ArgBiomass(idx_1),'Color',clr_ArgBiomass);

figure();
hold on;
plot(glc_list(idx_2),ref_mu(idx_2),'Color',clr_mu);
plot(glc_list(idx_2),ref_BiomassGlc(idx_2),'Color',clr_BiomassGlc);
plot(glc_list(idx_2),ref_ArgBiomass(idx_2),'Color',clr_ArgBiomass);

figure();
hold on;
plot(glc_list(idx_3),ref_mu(idx_3),'Color',clr_mu);
plot(glc_list(idx_3),ref_BiomassGlc(idx_3),'Color',clr_BiomassGlc);
plot(glc_list(idx_3),ref_ArgBiomass(idx_3),'Color',clr_ArgBiomass);
