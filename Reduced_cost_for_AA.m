%% Reduced cost for AA uptake
% Timing: ~ 34000 s

% Simulated results will be saved in the folder 'Results'.

tic;
load('pcLactis_Model.mat');
model = pcLactis_Model;

%% Optimization setting.

rxnID = 'R_dummy_assumed_Monomer';
osenseStr = 'Maximize';

%% Parameters.
GAM = 38;%ATP coefficient in the new biomass equation.
NGAM = 2; %(mmol/gCDW/h)
f_unmodeled = 0.45; %proportion of unmodeled protein in total protein (g/g)

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
model = changeRxnBounds(model,'R_M_PROTS_LLA_v3',0,'b');

% Block other glucose transporters
model = changeRxnBounds(model,'R_M_GLCpts_1',0,'b');
model = changeRxnBounds(model,'R_M_GLCt2_fwd',0,'b');

%% Main part.
% load glucose concentration of reference state
load('Sglc2_result_with_sf.mat');
selected_points = [1;5;9;13];
glc_conc_list = glc_conc_with_sf(selected_points,:);
clear glc_conc_with_sf;

AA_factor = 2;

% load AA id
aa_list = {'ala'
           'arg'
           'asn'
           'asp'
           'cys'
           'gln'
           'glu'
           'gly'
           'his'
           'ile'
           'leu'
           'lys'
           'met'
           'phe'
           'pro'
           'ser'
           'thr'
           'trp'
           'tyr'
           'val'};

% obtain the global saturation factor
load('Egsf2_result.mat');
x = global_saturation_factor_list(:,1);
y = global_saturation_factor_list(:,2);
y(3) = 1;
x = x(y ~= 1);
y = y(y ~= 1);
sf_coeff = x\y;
clear x y;

fluxes_rcAA = zeros(length(model.rxns),length(selected_points)*length(aa_list));

for i = 1:length(selected_points)
    
    mu_ref = glc_conc_list(i,1);
	glc_conc = glc_conc_list(i,2);
    factor_glc = glc_conc / (glc_conc + Km);

    for j = 1:length(aa_list)
        model_tmp = model;
        model_tmp = changeRxnBounds(model_tmp,Exchange_AAs,LBfactor_AAs*mu_ref,'l');
        
        aaid = aa_list(j);
        aarxnid = Exchange_AAs(contains(Exchange_AAs,aaid));
        aalb = LBfactor_AAs(contains(Exchange_AAs,aaid))*mu_ref*AA_factor;
        model_tmp = changeRxnBounds(model_tmp,aarxnid,aalb,'l');

        mu_low = 0;
        mu_high = 1;    
    
        while mu_high-mu_low > 0.001
            mu_mid = (mu_low+mu_high)/2;
            disp(['Glucose concentration = ' num2str(glc_conc) '; ' cell2mat(aaid) '; mu = ' num2str(mu_mid)]);
            model_tmptmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_mid,'b');
            
            factor_k = sf_coeff * mu_mid;
            if factor_k > 1
                factor_k = 1;
            end
            fileName = WriteLPSatFactor(model_tmptmp,mu_mid,f,osenseStr,rxnID,factor_k,...
                                        f_transporter,kcat_glc,factor_glc,...
                                        Info_enzyme,...
                                        Info_mRNA,...
                                        Info_protein,...
                                        Info_ribosome,...
                                        Info_tRNA);
            command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
            system(command,'-echo');
            fileName_out = 'Simulation.lp.out';
            [~,solME_status,~] = ReadSoplexResult(fileName_out,model_tmptmp);
            if strcmp(solME_status,'optimal')
                mu_low = mu_mid;
            else
                mu_high = mu_mid;
            end
        end

        model_tmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_low,'b');
        factor_k = sf_coeff * mu_low;
        if factor_k > 1
            factor_k = 1;
        end
        fileName = WriteLPSatFactor(model_tmp,mu_low,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,factor_glc,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,~,solME_full] = ReadSoplexResult(fileName_out,model_tmp);
        
        if strcmp(solME_status,'optimal')
            fluxes_rcAA(:,(i-1)*20+j) = solME_full;
        else
            fluxes_rcAA(:,(i-1)*20+j) = zeros(length(model_tmp.rxns),1);
        end
    end
end

% for i = 1:length(selected_points)
%     
%     mu_ref = glc_conc_list(i,1);
% 	glc_conc = glc_conc_list(i,2);
%     factor_glc = glc_conc / (glc_conc + Km);
% 
%     for j = 1:length(aa_list)
%         model_tmp = model;
%         aaid = aa_list(j);
%         
%         mu_low = 0;
%         mu_high = 2;    
%     
%         while mu_high-mu_low > 0.01
%             mu_mid = (mu_low+mu_high)/2;
%             disp(['Glucose concentration = ' num2str(glc_conc) cell2mat(aaid) '; mu = ' num2str(mu_mid)]);
%             model_tmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_mid,'b');
%             model_tmp = changeRxnBounds(model_tmp,Exchange_AAs,LBfactor_AAs*mu_mid,'l');
%             
%             aarxnid = Exchange_AAs(contains(Exchange_AAs,aaid));
%             aalb = LBfactor_AAs(contains(Exchange_AAs,aaid))*mu_mid*AA_factor;
%             model_tmp = changeRxnBounds(model_tmp,aarxnid,aalb,'l');
%             
%             factor_k = sf_coeff * mu_mid;
%             if factor_k > 1
%                 factor_k = 1;
%             end
%             fileName = WriteLPSatFactor(model_tmp,mu_mid,f,osenseStr,rxnID,factor_k,...
%                                         f_transporter,kcat_glc,factor_glc,...
%                                         Info_enzyme,...
%                                         Info_mRNA,...
%                                         Info_protein,...
%                                         Info_ribosome,...
%                                         Info_tRNA);
%             command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
%             system(command,'-echo');
%             fileName_out = 'Simulation.lp.out';
%             [~,solME_status,~] = ReadSoplexResult(fileName_out,model_tmp);
%             if strcmp(solME_status,'optimal')
%                 mu_low = mu_mid;
%             else
%                 mu_high = mu_mid;
%             end
%         end
%         
%         model_tmp = model;
%         model_tmp = changeRxnBounds(model_tmp,'R_biomass_dilution',mu_low,'b');
%         model_tmp = changeRxnBounds(model_tmp,Exchange_AAs,LBfactor_AAs*mu_low,'l');
%         aarxnid = Exchange_AAs(contains(Exchange_AAs,aaid));
%         aalb = LBfactor_AAs(contains(Exchange_AAs,aaid))*mu_low*AA_factor;
%         model_tmp = changeRxnBounds(model_tmp,aarxnid,aalb,'l');
%         
%         factor_k = sf_coeff * mu_low;
%         if factor_k > 1
%             factor_k = 1;
%         end
%         fileName = WriteLPSatFactor(model_tmp,mu_low,f,osenseStr,rxnID,factor_k,...
%                                     f_transporter,kcat_glc,factor_glc,...
%                                     Info_enzyme,...
%                                     Info_mRNA,...
%                                     Info_protein,...
%                                     Info_ribosome,...
%                                     Info_tRNA);
%         command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
%         system(command,'-echo');
%         fileName_out = 'Simulation.lp.out';
%         [~,~,solME_full] = ReadSoplexResult(fileName_out,model_tmp);
%         
%         if strcmp(solME_status,'optimal')
%             fluxes_rcAA(:,(i-1)*20+j) = solME_full;
%         else
%             fluxes_rcAA(:,(i-1)*20+j) = zeros(length(model_tmp.rxns),1);
%         end
%     end
% end

cd Results/;
save('RcAA_fluxes.mat','fluxes_rcAA');
cd ../;

clear;
toc;

%% Figures
load('RcAA_fluxes.mat');
flux_res = fluxes_rcAA;
clear fluxes_rcAA;

load('pcLactis_Model.mat');
model = pcLactis_Model;

mu = flux_res(strcmp(model.rxns,'R_biomass_dilution'),:);


