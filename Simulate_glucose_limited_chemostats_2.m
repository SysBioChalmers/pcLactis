%% Simulate glucose-limited chemostats (objective: minimizing glucose concentration)

% Timing: ~ 18600 s

% Without and with the saturation saturation factor are performed.

% Simulated results will be saved in the folder 'Results'.

% Figures can be obtained by running the codes starting from line 222.

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

% Block one of ADH isozymes llmg_0955
model = changeRxnBounds(model,'R_M_ALCD2x_1_rvs',0,'b');

%% Main part.

D_list = 0.1:0.05:0.7;%unit: /h

% without saturation factor
factor_k = 1;%global saturation factor
fluxes_simulated_without_sf = zeros(length(model.rxns),length(D_list));
glc_conc_without_sf = zeros(length(D_list),2);

for i = 1:length(D_list)
    
    D = D_list(i);
    
    factor_glc_low = 0;
    factor_glc_high = 1;
    
    while factor_glc_high-factor_glc_low > 0.001
        factor_glc_mid = (factor_glc_low+factor_glc_high)/2;
        disp(['Without sf: D = ' num2str(D) '; factor_glc = ' num2str(factor_glc_mid)]);
        model = changeRxnBounds(model,'R_biomass_dilution',D,'b');
        model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*D,'l');
        
        fileName = WriteLPSatFactor(model,D,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,factor_glc_mid,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);

        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,solME_status,~] = ReadSoplexResult(fileName_out,model);
        
        if strcmp(solME_status,'optimal')
            factor_glc_high = factor_glc_mid;
        else
            factor_glc_low = factor_glc_mid;
        end
    end
    
	fileName = WriteLPSatFactor(model,D,f,osenseStr,rxnID,factor_k,...
                                f_transporter,kcat_glc,factor_glc_high,...
                                Info_enzyme,...
                                Info_mRNA,...
                                Info_protein,...
                                Info_ribosome,...
                                Info_tRNA);
                   
	command =sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
	system(command,'-echo');
	fileName_out = 'Simulation.lp.out';
	[~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model);
    
    if strcmp(solME_status,'optimal')
        fluxes_simulated_without_sf(:,i) = solME_full;
        glc_conc = Km * factor_glc_high / (1 - factor_glc_high);
        glc_conc_without_sf(i,1) = D;
        glc_conc_without_sf(i,2) = glc_conc;
    else
        fluxes_simulated_without_sf(:,i) = zeros(length(model.rxns),1);
        glc_conc_without_sf(i,1) = D;
        glc_conc_without_sf(i,2) = glc_conc;
    end
end

% with saturation factor
% obtain the global saturation factor
load('Egsf2_result.mat');
x = global_saturation_factor_list(:,1);
y = global_saturation_factor_list(:,2);
% y(3) = 1;
x = x(y ~= 1);
y = y(y ~= 1);
sf_coeff = x\y;
clear x y;

fluxes_simulated_with_sf = zeros(length(model.rxns),length(D_list));
glc_conc_with_sf = zeros(length(D_list),2);

for i = 1:length(D_list)
    
    D = D_list(i);
    
	factor_k = sf_coeff * D;
    if factor_k > 1
        factor_k = 1;
    end
    
    factor_glc_low = 0;
    factor_glc_high = 1;
    
    while factor_glc_high-factor_glc_low > 0.001
        factor_glc_mid = (factor_glc_low+factor_glc_high)/2;
        disp(['With sf: D = ' num2str(D) '; factor_glc = ' num2str(factor_glc_mid)]);
        model = changeRxnBounds(model,'R_biomass_dilution',D,'b');
        model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*D,'l');
        
        fileName = WriteLPSatFactor(model,D,f,osenseStr,rxnID,factor_k,...
                                    f_transporter,kcat_glc,factor_glc_mid,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);

        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,solME_status,~] = ReadSoplexResult(fileName_out,model);
        
        if strcmp(solME_status,'optimal')
            factor_glc_high = factor_glc_mid;
        else
            factor_glc_low = factor_glc_mid;
        end
    end
    
	fileName = WriteLPSatFactor(model,D,f,osenseStr,rxnID,factor_k,...
                                f_transporter,kcat_glc,factor_glc_high,...
                                Info_enzyme,...
                                Info_mRNA,...
                                Info_protein,...
                                Info_ribosome,...
                                Info_tRNA);
                   
	command =sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
	system(command,'-echo');
	fileName_out = 'Simulation.lp.out';
	[~,solME_status,solME_full] = ReadSoplexResult(fileName_out,model);
    
    if strcmp(solME_status,'optimal')
        fluxes_simulated_with_sf(:,i) = solME_full;
        glc_conc = Km * factor_glc_high / (1 - factor_glc_high);
        glc_conc_with_sf(i,1) = D;
        glc_conc_with_sf(i,2) = glc_conc;
    else
        fluxes_simulated_with_sf(:,i) = zeros(length(model.rxns),1);
        glc_conc_with_sf(i,1) = D;
        glc_conc_with_sf(i,2) = 0;
    end
end

cd Results/;
save('Sglc2_fluxes_without_sf.mat','fluxes_simulated_without_sf');
save('Sglc2_fluxes_with_sf.mat','fluxes_simulated_with_sf');
save('Sglc2_result_without_sf.mat','glc_conc_without_sf');
save('Sglc2_result_with_sf.mat','glc_conc_with_sf');
cd ../;

clear;
toc;

%% Figure.
load('Sglc2_fluxes_without_sf.mat');
load('Sglc2_fluxes_with_sf.mat');
load('Sglc2_result_without_sf.mat');
load('Sglc2_result_with_sf.mat');

% %Glucose concentration and growth rate
% figure('Name','1');
% hold on;
% box on;
% y1 = glc_conc_without_sf(:,1);
% x1 = glc_conc_without_sf(:,2);
% plot(x1,y1,'-o','LineWidth',1,'Color',[69,117,180]/255,'MarkerSize',8,'MarkerEdgeColor',[69,117,180]/255);
% y2 = glc_conc_with_sf(:,1);
% x2 = glc_conc_with_sf(:,2);
% plot(x2,y2,'-o','LineWidth',1,'Color',[215,48,39]/255,'MarkerSize',8,'MarkerEdgeColor',[215,48,39]/255);
% legend({'Before','After'},'FontSize',12,'FontName','Helvetica','location','north');
% set(gca,'FontSize',12,'FontName','Helvetica');
% xlabel('Glucose concentration (uM)','FontSize',14,'FontName','Helvetica');
% ylabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
% 
% set(gcf,'position',[0 300 230 130]);
% set(gca,'position',[0.23 0.25 0.72 0.7]);
% clear x1 x2 y1 y2;

%Metabolic shift
% Experimental data
[~, ~, exp_raw] = xlsread('Exchange_reaction_setting.xlsx','Exchange_rates');
exp_mu = cell2mat(exp_raw(2,14:18));

f_exp_mix = cell2mat(exp_raw(15,14:18));
f_exp_lac = cell2mat(exp_raw(14,14:18));
f_exp_mix_sd = cell2mat(exp_raw(15,19:23));
f_exp_lac_sd = cell2mat(exp_raw(14,19:23));

ave_orn = cell2mat(exp_raw(8,14:18));
ave_arg = cell2mat(exp_raw(9,14:18));
ave_nh4 = cell2mat(exp_raw(10,14:18));
ave_citr = cell2mat(exp_raw(11,14:18));
sd_orn = cell2mat(exp_raw(8,19:23));
sd_arg = cell2mat(exp_raw(9,19:23));
sd_nh4 = cell2mat(exp_raw(10,19:23));
sd_citr = cell2mat(exp_raw(11,19:23));

% Simulated data
load('pcLactis_Model.mat');
model = pcLactis_Model;
mu1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_biomass_dilution'),:);
glc1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_glc__D_e'),:);
ac1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_ac_e'),:);
eth1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_etoh_e'),:);
form1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_for_e'),:);
lac1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_lac__L_e'),:);
orn1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_orn_e'),:);
citr1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_citr__L_e'),:);
nh41 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_nh4_e'),:);
arg1 = fluxes_simulated_without_sf(strcmp(model.rxns,'R_M_EX_arg__L_e'),:);
f_mix1 = (ac1*2+eth1*2+form1)./(-glc1*6);
f_lac1 = (lac1*3)./(-glc1*6);

mu2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_biomass_dilution'),:);
glc2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_glc__D_e'),:);
ac2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_ac_e'),:);
eth2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_etoh_e'),:);
form2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_for_e'),:);
lac2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_lac__L_e'),:);
orn2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_orn_e'),:);
citr2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_citr__L_e'),:);
nh42 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_nh4_e'),:);
arg2 = fluxes_simulated_with_sf(strcmp(model.rxns,'R_M_EX_arg__L_e'),:);
f_mix2 = (ac2*2+eth2*2+form2)./(-glc2*6);
f_lac2 = (lac2*3)./(-glc2*6);

% Figures
% mixed to lactic acid without sf
color_mixed = [94,60,153]/255;
color_lactic = [230,97,1]/255;
figure('Name','2');
hold on;
box on;
errorbar(exp_mu,f_exp_lac,f_exp_lac_sd,'-.o','LineWidth',0.75,'Color',color_lactic,'MarkerSize',8,'MarkerEdgeColor',color_lactic);
errorbar(exp_mu,f_exp_mix,f_exp_mix_sd,'-.o','LineWidth',0.75,'Color',color_mixed,'MarkerSize',8,'MarkerEdgeColor',color_mixed);
plot(mu1,f_lac1,'-o','LineWidth',0.75,'Color',color_lactic,'MarkerSize',8,'MarkerEdgeColor',color_lactic,'MarkerFaceColor',color_lactic);
plot(mu1,f_mix1,'-o','LineWidth',0.75,'Color',color_mixed,'MarkerSize',8,'MarkerEdgeColor',color_mixed,'MarkerFaceColor',color_mixed);
set(gca,'YTick',0:0.2:1);
set(gca,'ycolor','k');
set(gca,'XTick',0.1:0.2:0.7);
xlim([0.1 0.7]);
ylim([-0.05 1.05]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Fraction in total glucose carbon','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
% legend({'Sim mixed acid',...
%         'Sim lactate'...
%         'Exp mixed acid',...
%         'Exp lactate'},'FontSize',12,'FontName','Helvetica','location','north');
set(gcf,'position',[0 0 305 300]);
set(gca,'position',[0.17 0.13 0.7 0.7]);

% mixed to lactic acid with sf
figure('Name','3');
hold on;
box on;

plot(mu2,f_lac2,'-','LineWidth',0.75,'Color',color_lactic);
x = exp_mu; y = f_exp_lac; yu = f_exp_lac + f_exp_lac_sd; yl = f_exp_lac - f_exp_lac_sd;
plot(x,y,'-.','LineWidth',0.75,'Color',color_lactic);
fill([x fliplr(x)],[yu fliplr(yl)],color_lactic,'linestyle','none','FaceAlpha',0.3);

plot(mu2,f_mix2,'-','LineWidth',0.75,'Color',color_mixed);
x = exp_mu; y = f_exp_mix; yu = f_exp_mix + f_exp_mix_sd; yl = f_exp_mix - f_exp_mix_sd;
plot(x,y,'-.','LineWidth',0.75,'Color',color_mixed);
fill([x fliplr(x)],[yu fliplr(yl)],color_mixed,'linestyle','none','FaceAlpha',0.3);

% errorbar(exp_mu,f_exp_lac,f_exp_lac_sd,'-.','LineWidth',0.75,'Color',color_lactic);
% errorbar(exp_mu,f_exp_mix,f_exp_mix_sd,'-.','LineWidth',0.75,'Color',color_mixed);

set(gca,'YTick',0:0.2:1);
set(gca,'ycolor','k');
set(gca,'XTick',0.1:0.2:0.7);
xlim([0.1 0.7]);
ylim([-0.05 1.05]);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Fraction in glucose','FontSize',14,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
legend({'Sim lactic acid',...
        'Exp lactic acid'...
        'SD lactic acid'...
        'Sim mixed acid',...
        'Exp mixed acid'...
        'SD mixed acid'},'FontSize',12,'FontName','Helvetica','location','north');
set(gcf,'position',[450 0 305 300]);
set(gca,'position',[0.17 0.13 0.7 0.7]);

% mixed to lactic acid with sf
figure('Name','3');
hold on;
box on;

plot(mu2,f_lac2,'-','LineWidth',0.75,'Color',color_lactic);
plot(mu2,f_mix2,'-','LineWidth',0.75,'Color',color_mixed);
set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Simulated growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel('Fraction in total glucose carbon','FontSize',14,'FontName','Helvetica');
set(gcf,'position',[450 0 305 300]);
set(gca,'position',[0.17 0.13 0.7 0.7]);

ax1 = gca;
ax1.XLim = [0.2 0.7];
set(ax1,'XTick',0.2:0.1:0.7);
ax1.YLim = [0 1];
ax2 = axes('Position',ax1.Position,'XAxisLocation','top','YAxisLocation','right','Color','none');
hold on;

x = exp_mu; y = f_exp_lac; yu = f_exp_lac + f_exp_lac_sd; yl = f_exp_lac - f_exp_lac_sd;
plot(x,y,'-.','LineWidth',0.75,'Color',color_lactic);
fill([x fliplr(x)],[yu fliplr(yl)],color_lactic,'linestyle','none','FaceAlpha',0.3);

x = exp_mu; y = f_exp_mix; yu = f_exp_mix + f_exp_mix_sd; yl = f_exp_mix - f_exp_mix_sd;
plot(x,y,'-.','LineWidth',0.75,'Color',color_mixed);
fill([x fliplr(x)],[yu fliplr(yl)],color_mixed,'linestyle','none','FaceAlpha',0.3);

set(gca,'XAxisLocation','top','YAxisLocation','right','Color','none');
set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Experimental growth rate (/h)','FontSize',14,'FontName','Helvetica');
ax2.XLim = [0.15 0.6];
ax2.YLim = ax1.YLim;
line(xlim,[ax2.YLim(2), ax2.YLim(2)],'color','w','LineWidth',1);
line(xlim,[ax2.YLim(2), ax2.YLim(2)],'color','k','LineStyle','-.');
legend({'lactic acid','lactic acid','mixed acid','mixed acid'},'FontSize',12,'FontName','Helvetica','location','nw');
set(ax2,'XTick',0.2:0.1:0.6);
set(ax2,'ytick',[]);

% arg catabolism with sf
color_orn = [221,52,151]/255;
color_arg = [65,171,93]/255;
color_citr = [82,82,82]/255;
figure('Name','4');
hold on;
plot(mu2,orn2,'-','LineWidth',0.75,'Color',color_orn);
plot(mu2,arg2,'-','LineWidth',0.75,'Color',color_arg);
plot(mu2,citr2,'-','LineWidth',0.75,'Color',color_citr);
set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Simulated growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel('Flux (mmol/gCDW/h)','FontSize',14,'FontName','Helvetica');
set(gcf,'position',[900 0 305 300]);
set(gca,'position',[0.17 0.13 0.7 0.7]);
ax1 = gca;
ax1.XLim = [0.2 0.7];
set(ax1,'XTick',0.2:0.1:0.7);
ax1.YLim = [-1.1 1.1];
ax2 = axes('Position',ax1.Position,'XAxisLocation','top','YAxisLocation','right','Color','none');
hold on;

x = exp_mu; y = ave_orn; yu = ave_orn + sd_orn; yl = ave_orn - sd_orn;
plot(x,y,'-.','LineWidth',0.75,'Color',color_orn);
fill([x fliplr(x)],[yu fliplr(yl)],color_orn,'linestyle','none','FaceAlpha',0.3);

x = exp_mu; y = ave_arg; yu = ave_arg + sd_arg; yl = ave_arg - sd_arg;
plot(x,y,'-.','LineWidth',0.75,'Color',color_arg);
fill([x fliplr(x)],[yu fliplr(yl)],color_arg,'linestyle','none','FaceAlpha',0.3);

x = exp_mu; y = ave_citr; yu = ave_citr + sd_citr; yl = ave_citr - sd_citr;
plot(x,y,'-.','LineWidth',0.75,'Color',color_citr);
fill([x fliplr(x)],[yu fliplr(yl)],color_citr,'linestyle','none','FaceAlpha',0.3);

% errorbar(exp_mu,ave_orn,sd_orn,'-.o','LineWidth',0.75,'Color',color_orn,'MarkerSize',8,'MarkerEdgeColor',color_orn);
% errorbar(exp_mu,ave_arg,sd_arg,'-.o','LineWidth',0.75,'Color',color_arg,'MarkerSize',8,'MarkerEdgeColor',color_arg);
% errorbar(exp_mu,ave_citr,sd_citr,'-.o','LineWidth',0.75,'Color',color_citr,'MarkerSize',8,'MarkerEdgeColor',color_citr);

set(gca,'XAxisLocation','top','YAxisLocation','right','Color','none');
set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Experimental growth rate (/h)','FontSize',14,'FontName','Helvetica');
ax2.XLim = [0.15 0.6];
ax2.YLim = ax1.YLim;
line(xlim,[ax2.YLim(2), ax2.YLim(2)],'color','w','LineWidth',1);
line(xlim,[ax2.YLim(2), ax2.YLim(2)],'color','k','LineStyle','-.');
legend({'Ornithine','Ornithine','Arginine','Arginine','Citrulline','Citrulline'},'FontSize',12,'FontName','Helvetica','location','nw');
set(ax2,'XTick',0.2:0.1:0.6);
set(ax2,'ytick',[]);

clear;
