%% Estimate GAM and NGAM for pcModel

% Timing: ~ 450 s

% We used physiological data under chemostat conditions (PMID: 25828364) as
% constraints for exchange reactions to estimate max ATP production rates.

tic;
load('pcLactis_Model.mat');
model = pcLactis_Model;

%% Optimization setting.
rxnID = 'R_M_ATPM';
osenseStr = 'Maximize';

%% Parameters.
GAM = 0;%ATP coefficient in the new biomass equation.
f_unmodeled = 0.4; %proportion of unmodeled protein in total protein (g/g)

model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',0,'l');
model = changeRxnBounds(model,'R_M_ATPM',1000,'u');
[model,f] = ChangeUnmodeledProtein(model,f_unmodeled);

kcat_glc = 180;%kcat value of glucose transporter
factor_k = 1;%global saturation factor
f_transporter = 1;%fraction of glucose transporter in total proteome

%% Data import.
load('Info_enzyme.mat');
load('Info_mRNA.mat');
load('Info_protein.mat');
load('Info_ribosome.mat');
load('Info_tRNA.mat');

%% Set reaction rates.
% Block uptake direction of all the exchange reactions.
idx = cellfun(@(x) contains(x,'R_M_EX_'),model.rxns,'UniformOutput',false);
M_exchange_reactions = model.rxns(cell2mat(idx));
model = changeRxnBounds(model,M_exchange_reactions,0,'l');
clear idx M_exchange_reactions;

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

%% Loop for dilution rate of 0.15 0.3 0.45 0.5 and 0.6.
[~, ~, exchange_raw] = xlsread('Exchange_reaction_setting.xlsx','Exp_bounds_GAM_NGAM');

Exchange_reactions = exchange_raw(2:end,1);
header = exchange_raw(1,1:21);

mu_list = [0.15 0.3 0.45 0.5 0.6];

res = zeros(0,2); %mu and q_atp list
fluxes = zeros(length(model.rxns),0);
for i = 1:length(mu_list)
    mu = mu_list(i);
    
    model_tmp = changeRxnBounds(model,'R_biomass_dilution',mu,'b');

    % Set bounds for some exchange reactions.
    idx = find(contains(header,num2str(mu)));
    replicate = length(idx)/2;
    
	for j = 1:replicate
        lb_idx = idx(2*j-1);
        ub_idx = idx(2*j);
        LB = cell2mat(exchange_raw(2:end,lb_idx));
        UB = cell2mat(exchange_raw(2:end,ub_idx));
        model_tmp = changeRxnBounds(model_tmp,Exchange_reactions,LB,'l');
        model_tmp = changeRxnBounds(model_tmp,Exchange_reactions,UB,'u');

        disp(['mu = ' num2str(mu) '; replicate = ' num2str(j)]);        
        
        fileName = WriteLPSatFactor(model_tmp,mu,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,1,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-15 -o1e-15 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [solME_obj,solME_status,solME_full] = ReadSoplexResult(fileName_out,model_tmp);
        
        if strcmp(solME_status,'optimal')
            fluxes = cat(2,fluxes,solME_full);
            res = [res;mu solME_obj];
        else
            fluxes = cat(2,fluxes,zeros(length(model_tmp.rxns),1));
        end
	end
end

x = res(:,1);
y = res(:,2);
p = polyfit(x,y,1);
GAM = p(1);
NGAM = p(2);
yfit =  p(1) * x + p(2);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;

clear x y p yfit yresid SSresid SStotal;

clear ans command exchange_raw f f_transporter f_unmodeled factor_k Exchange_reactions;
clear fileName fileName_out GAM header i idx j kcat_glc LB lb_idx pcLactis_Model;
clear model_tmp mu mu_list NGAM osenseStr replicate rxnID solME_status UB ub_idx;
clear Info_enzyme Info_mRNA Info_protein Info_protein Info_ribosome Info_tRNA solME_full solME_obj;

toc;