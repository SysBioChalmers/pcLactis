%% Sensitivity analysis for total modeled proteome and glucose transporter.

% Timing: ~ 60000 s

% With the saturation saturation factor are performed.

% Simulated results will be saved in the folder 'Results'.

% Figures can be obtained by running the codes starting from line .

% tic;
% load('pcLactis_Model.mat');
% model = pcLactis_Model;
% 
% %% Optimization setting.
% 
% rxnID = 'R_dummy_assumed_Monomer';
% osenseStr = 'Maximize';
% 
% %% Parameters.
% GAM = 42;%ATP coefficient in the new biomass equation.
% NGAM = 2.5; %(mmol/gCDW/h)
% 
% model = ChangeATPinBiomass(model,GAM);
% model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
% 
% kcat_glc = 180;%kcat value of glucose transporter
% Km = 21;%Km of glucose transporter, unit: uM (PMID: 30630406)
% 
% %% Data import.
% load('Info_enzyme.mat');
% load('Info_mRNA.mat');
% load('Info_protein.mat');
% load('Info_ribosome.mat');
% load('Info_tRNA.mat');
% [~, ~, exchange_raw] = xlsread('Exchange_reaction_setting.xlsx','MaxMu');
% [~, ~, aa_raw] = xlsread('Exchange_reaction_setting.xlsx','AA_factors');
% 
% Exchange_AAs = aa_raw(2:end,1);
% LBfactor_AAs = cell2mat(aa_raw(2:end,2));
% clear aa_raw;
% 
% %% Set reaction rates.
% % Block uptake direction of all the exchange reactions.
% idx = cellfun(@(x) contains(x,'R_M_EX_'),model.rxns,'UniformOutput',false);
% M_exchange_reactions = model.rxns(cell2mat(idx));
% model = changeRxnBounds(model,M_exchange_reactions,0,'l');
% clear idx M_exchange_reactions;
% % Set bounds for some exchange reactions.
% Exchange_reactions = exchange_raw(2:end,1);
% LB = cell2mat(exchange_raw(2:end,2));
% UB = cell2mat(exchange_raw(2:end,3));
% model = changeRxnBounds(model,Exchange_reactions,LB,'l');
% model = changeRxnBounds(model,Exchange_reactions,UB,'u');
% clear exchange_raw Exchange_reactions LB UB;
% % Block some reactions in the M model.
% model = changeRxnBounds(model,'R_M_biomass_LLA',0,'b');
% model = changeRxnBounds(model,'R_M_biomass_LLA_atpm',0,'b');
% model = changeRxnBounds(model,'R_M_PROTS_LLA',0,'b');
% model = changeRxnBounds(model,'R_M_PROTS_LLA_v2',0,'b');
% model = changeRxnBounds(model,'R_M_PROTS_LLA_v3',0,'b');
% model = changeRxnBounds(model,'R_M_MGt2pp_rvs',0,'b');%block infinite h[e]
% 
% % Block other glucose transporters
% model = changeRxnBounds(model,'R_M_GLCpts_2',0,'b');
% model = changeRxnBounds(model,'R_M_GLCpermease_fwd',0,'b');
% 
% %% Main part.
% 
% glc_list = [2 4 6 8 10 20 40 60 80 100 200 1000 10000 100000 1000000];%unit: /uM
% modeled_protein_list = (1-0.42)*(1:0.1:1.5);
% glucose_transporter_list = 0.01*(1:0.1:1.5);
% res1 = zeros(length(modeled_protein_list),4,length(glc_list));
% res2 = zeros(length(glucose_transporter_list),4,length(glc_list));
% 
% % obtain the global saturation factor
% load('Egsf1_result.mat');
% x = global_saturation_factor_list(:,1);
% y = global_saturation_factor_list(:,2);
% x = x(y ~= 1);
% y = y(y ~= 1);
% sf_coeff = x\y;
% clear x y;
% 
% for i = 1:length(glc_list)
%     glc_conc = glc_list(i);
%     factor_glc = glc_conc / (glc_conc + Km);
%     
%     for j = 1:length(modeled_protein_list)
%         f_unmodeled = 1-modeled_protein_list(j);
%         [model,f] = ChangeUnmodeledProtein(model,f_unmodeled);
%         f_transporter = 0.01;
%         
%         mu_low = 0;
%         mu_high = 2;
% 
%         while mu_high-mu_low > 0.01
%             mu_mid = (mu_low+mu_high)/2;
%             disp(['Glucose concentration = ' num2str(glc_conc) ' uM; modeled protein = ' num2str((1-f_unmodeled)) '; glucose transporter = ' num2str(f_transporter) '; mu = ' num2str(mu_mid)]);
%             model = changeRxnBounds(model,'R_biomass_dilution',mu_mid,'b');
%             model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*mu_mid,'l');
%             factor_k = sf_coeff * mu_mid;
%                 if factor_k > 1
%                     factor_k = 1;
%                 end
%             fileName = WriteLPSatFactor(model,mu_mid,f,osenseStr,rxnID,factor_k,...
%                                         f_transporter,kcat_glc,factor_glc,...
%                                         Info_enzyme,...
%                                         Info_mRNA,...
%                                         Info_protein,...
%                                         Info_ribosome,...
%                                         Info_tRNA);
%             command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
%             system(command,'-echo');
%             fileName_out = 'Simulation.lp.out';
%             [~,solME_status,~] = ReadSoplexResult(fileName_out,model);
%             if strcmp(solME_status,'optimal')
%                 mu_low = mu_mid;
%             else
%                 mu_high = mu_mid;
%             end
%         end
%         
%         model = changeRxnBounds(model,'R_biomass_dilution',mu_low,'b');
%         model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*mu_low,'l');
%         factor_k = sf_coeff * mu_low;
%         if factor_k > 1
%             factor_k = 1;
%         end
%         fileName = WriteLPSatFactor(model,mu_low,f,osenseStr,rxnID,factor_k,...
%                                     f_transporter,kcat_glc,factor_glc,...
%                                     Info_enzyme,...
%                                     Info_mRNA,...
%                                     Info_protein,...
%                                     Info_ribosome,...
%                                     Info_tRNA);
%         command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
%         system(command,'-echo');
%         fileName_out = 'Simulation.lp.out';
%         [~,~,solME_full] = ReadSoplexResult(fileName_out,model);
%         
%         mu = solME_full(strcmp(model.rxns,'R_biomass_dilution'),1);
%         glc = -solME_full(strcmp(model.rxns,'R_M_EX_glc_LPAREN_e_RPAREN_'),1);
%         arg = -solME_full(strcmp(model.rxns,'R_M_EX_arg_L_LPAREN_e_RPAREN_'),1);
%         res1(j,:,i) = [1-f_unmodeled mu mu/(glc*180/1000) arg*174/1000/mu];
%                                        %g_CDW/g_glucose   g_arginine/g_CDW
%     end
%         
% 	for k = 1:length(glucose_transporter_list)
%         f_transporter = glucose_transporter_list(k);
%         f_unmodeled = 0.42;
%         [model,f] = ChangeUnmodeledProtein(model,f_unmodeled);
%         
%         mu_low = 0;
%         mu_high = 2;
%         
%         while mu_high-mu_low > 0.01
%             mu_mid = (mu_low+mu_high)/2;
%             disp(['Glucose concentration = ' num2str(glc_conc) ' uM; modeled protein = ' num2str((1-f_unmodeled)) '; glucose transporter = ' num2str(f_transporter) '; mu = ' num2str(mu_mid)]);
%             model = changeRxnBounds(model,'R_biomass_dilution',mu_mid,'b');
%             model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*mu_mid,'l');
%             factor_k = sf_coeff * mu_mid;
%             if factor_k > 1
%                 factor_k = 1;
%             end
%             fileName = WriteLPSatFactor(model,mu_mid,f,osenseStr,rxnID,factor_k,...
%                                         f_transporter,kcat_glc,factor_glc,...
%                                         Info_enzyme,...
%                                         Info_mRNA,...
%                                         Info_protein,...
%                                         Info_ribosome,...
%                                         Info_tRNA);
%             command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
%             system(command,'-echo');
%             fileName_out = 'Simulation.lp.out';
%             [~,solME_status,~] = ReadSoplexResult(fileName_out,model);
%             if strcmp(solME_status,'optimal')
%                 mu_low = mu_mid;
%             else
%                 mu_high = mu_mid;
%             end
%         end
%         
%         model = changeRxnBounds(model,'R_biomass_dilution',mu_low,'b');
%         model = changeRxnBounds(model,Exchange_AAs,LBfactor_AAs*mu_low,'l');
%         factor_k = sf_coeff * mu_low;
%         if factor_k > 1
%             factor_k = 1;
%         end
%         fileName = WriteLPSatFactor(model,mu_low,f,osenseStr,rxnID,factor_k,...
%                                     f_transporter,kcat_glc,factor_glc,...
%                                     Info_enzyme,...
%                                     Info_mRNA,...
%                                     Info_protein,...
%                                     Info_ribosome,...
%                                     Info_tRNA);
%         command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-10 -o1e-10 -x -q -c --readmode=1 --solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
%         system(command,'-echo');
%         fileName_out = 'Simulation.lp.out';
%         [~,~,solME_full] = ReadSoplexResult(fileName_out,model);
%         
%         mu = solME_full(strcmp(model.rxns,'R_biomass_dilution'),1);
%         glc = -solME_full(strcmp(model.rxns,'R_M_EX_glc_LPAREN_e_RPAREN_'),1);
%         arg = -solME_full(strcmp(model.rxns,'R_M_EX_arg_L_LPAREN_e_RPAREN_'),1);
%         res2(k,:,i) = [f_transporter mu mu/(glc*180/1000) arg*174/1000/mu];
%                                        %g_CDW/g_glucose   g_arginine/g_CDW
% 	end
% end
%     
% 
% 
% cd Results/;
% save('Ac_result1.mat','res1');
% save('Ac_result2.mat','res2');
% cd ../;
% 
% clear;
% toc;

%% Data import
load('Ac_result1.mat');
load('Ac_result2.mat');
glc_list = [2 4 6 8 10 20 40 60 80 100 200 1000 10000 100000 1000000];%unit: /uM
idx_1 = 2:4;
idx_2 = 5:10;
idx_3 = 11:length(glc_list);
st1 = glc_list(idx_1);%unit: /uM
st2 = glc_list(idx_2);%unit: /uM
st3 = glc_list(idx_3);%unit: /uM

clr_mu = [227,26,28]/255;
clr_BiomassGlc = [31,120,180]/255;
clr_ArgBiomass = [51,160,44]/255;

%% Figure for reference
ref_mu = permute(res1(1,2,:),[3 2 1]);
ref_BiomassGlc = permute(res1(1,3,:),[3 2 1]);
ref_ArgBiomass = permute(res1(1,4,:),[3 2 1]);

% stage 1
figure('position',[0 0 800 300]);
pos_st1 = [0.2 0.13 0.06 0.7];
h1 = axes('position',pos_st1,'color','none','ycolor',clr_mu);
hold on;
plot(glc_list(idx_1),ref_mu(idx_1),'Color',clr_mu);
set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Glucose concentration (uM)','FontSize',14,'FontName','Helvetica');
ylabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
h1.YLim = [0 0.8];
set(h1,'ytick',0:0.2:0.8);
h1.XLim = [min(glc_list(idx_1)) max(glc_list(idx_1))];
h2 = axes('Position',h1.Position,'Color','none');
hold on;
plot(glc_list(idx_1),ref_BiomassGlc(idx_1),'Color',clr_BiomassGlc);
h2.YLim = [0 0.3];
h2.XLim = [min(glc_list(idx_1)) max(glc_list(idx_1))];
set(h2,'ycolor',clr_mu,'xtick',[],'ytick',[]);
h3 = axes('position',[0.14 0.13 0.01 0.7],'color','none','ycolor',clr_BiomassGlc);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Biomass yield on glucose (g/g)','FontSize',14,'FontName','Helvetica');
set(h3,'yaxislocation','left','xtick',[]);
h3.YLim = h2.YLim;
set(h3,'ytick',0:0.1:0.3);
h4 = axes('Position',h1.Position,'Color','none');
hold on;
plot(glc_list(idx_1),ref_ArgBiomass(idx_1),'Color',clr_ArgBiomass);
h4.YLim = [0 0.3];
h4.XLim = [min(glc_list(idx_1)) max(glc_list(idx_1))];
set(h4,'ycolor',clr_mu,'xtick',[],'ytick',[]);
h5 = axes('position',[0.08 0.13 0.01 0.7],'color','none','ycolor',clr_ArgBiomass);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Arg consumed per biomass (g/g)','FontSize',14,'FontName','Helvetica');
set(h5,'yaxislocation','left','xtick',[]);
h5.YLim = h4.YLim;
set(h5,'ytick',0:0.1:0.3);

% stage 2
pos_st2 = [0.29 0.13 0.11 0.7];
h6 = axes('position',pos_st2,'color','none');
hold on;
plot(glc_list(idx_2),ref_mu(idx_2),'Color',clr_mu);
h6.YLim = [0 0.8];
h6.XLim = [min(glc_list(idx_2)) max(glc_list(idx_2))];
axis off;
h7 = axes('position',[0.29 0.13 0.11 0.01],'color','none');
set(gca,'FontSize',12,'FontName','Helvetica');
set(h7,'ytick',[]);
h7.XLim = h6.XLim;
set(h7,'xtick',[min(glc_list(idx_2)),max(glc_list(idx_2))/2,max(glc_list(idx_2))]);
h8 = axes('position',pos_st2,'color','none');
hold on;
plot(glc_list(idx_2),ref_BiomassGlc(idx_2),'Color',clr_BiomassGlc);
h8.YLim = [0 0.3];
h8.XLim = h6.XLim;
axis off;
h9 = axes('position',pos_st2,'color','none');
hold on;
plot(glc_list(idx_2),ref_ArgBiomass(idx_2),'Color',clr_ArgBiomass);
h9.YLim = [0 0.3];
h9.XLim = h6.XLim;
axis off;

% stage 3
pos_st3 = [0.43 0.13 0.08 0.7];
h10 = axes('position',pos_st3,'color','none');
set(gca,'xscale','log');
hold on;
plot(glc_list(idx_3),ref_mu(idx_3),'Color',clr_mu);
h10.YLim = [0 0.8];
h10.XLim = [min(glc_list(idx_3)) max(glc_list(idx_3))];
axis off;
h11 = axes('position',[0.43 0.13 0.08 0.01],'color','none');
set(gca,'FontSize',12,'FontName','Helvetica','xscale','log');
set(h11,'ytick',[]);
h11.XLim = h10.XLim;
h12 = axes('position',pos_st3,'color','none');
set(gca,'xscale','log');
hold on;
plot(glc_list(idx_3),ref_BiomassGlc(idx_3),'Color',clr_BiomassGlc);
h12.YLim = [0 0.3];
h12.XLim = h10.XLim;
axis off;
h13 = axes('position',pos_st3,'color','none');
hold on;
set(gca,'xscale','log');
plot(glc_list(idx_3),ref_ArgBiomass(idx_3),'Color',clr_ArgBiomass);
h13.YLim = [0 0.3];
h13.XLim = h10.XLim;
axis off;

%% Figures for analysing growth rates
figure('position',[0 0 1200 600]);
az = -20;
el = 20;

% relative data for changing modeled protein
[a,b,c] = size(res1);
rel_res1 = zeros(a,b,c);
for i = 1:c
    tmp = res1(:,:,i);
    rel_res1(:,:,i) = tmp./tmp(1,:);
end
% relative data for changing glucose transporter
[a,b,c] = size(res2);
rel_res2 = zeros(a,b,c);
for i = 1:c
    tmp = res2(:,:,i);
    rel_res2(:,:,i) = tmp./tmp(1,:);
end

idx = [idx_1 idx_2 idx_3];
max_tmp1 = max(rel_res1(:,2,:));
max_tmp1 = max_tmp1(idx);
max_tmp2 = max(rel_res2(:,2,:));
max_tmp2 = max_tmp2(idx);
min_tmp1 = min(rel_res1(:,2,:));
min_tmp1 = min_tmp1(idx);
min_tmp2 = min(rel_res2(:,2,:));
min_tmp2 = min_tmp2(idx);
max_rel_mu = max([max(max_tmp1),max(max_tmp2)]);
min_rel_mu = min([min(min_tmp1),min(min_tmp2)]);
range_mu = max_rel_mu-min_rel_mu;

clrmax = clr_mu;
n_points = 64;

% changing modeled protein
% stage 1
ax1 = subplot(2,3,1);
data1 = rel_res1(:,:,idx_1);
x = st1.*ones(a,length(idx_1));
y = permute(data1(:,1,:),[1 3 2]);
z1 = permute(data1(:,2,:),[1 3 2]);
surf(x,y,z1,'EdgeColor',[1 1 1],'FaceColor','interp');
xlim([min(st1),max(st1)]);
ylim([1,1.5]);
zlim([1,1.7]);
set(ax1,'xtick',[min(st1),max(st1)]);
set(ax1,'ztick',[1,1.5]);
view(az, el);
max_rel_mu_tmp = max(max(data1(:,2,:)));
min_rel_mu_tmp = min(min(data1(:,2,:)));
max_clr = (max_rel_mu_tmp-min_rel_mu)/range_mu*clrmax;
min_clr = (min_rel_mu_tmp-min_rel_mu)/range_mu*clrmax;
clr = [linspace(min_clr(1),max_clr(1),n_points);linspace(min_clr(2),max_clr(2),n_points);linspace(min_clr(3),max_clr(3),n_points)]';
colormap(ax1,clr);
colorbar;
set(ax1,'FontSize',12,'FontName','Helvetica');

% stage 2
ax2 = subplot(2,3,2);
data2 = rel_res1(:,:,idx_2);
x = st2.*ones(a,length(idx_2));
y = permute(data2(:,1,:),[1 3 2]);
z1 = permute(data2(:,2,:),[1 3 2]);
surf(x,y,z1,'EdgeColor',[1 1 1],'FaceColor','interp');
xlim([min(st2),max(st2)]);
ylim([1,1.5]);
zlim([1,1.7]);
set(ax2,'xtick',[min(st2),max(st2)]);
set(ax2,'ztick',[1,1.5]);
view(az, el);
max_rel_mu_tmp = max(max(data2(:,2,:)));
min_rel_mu_tmp = min(min(data2(:,2,:)));
max_clr = (max_rel_mu_tmp-min_rel_mu)/range_mu*clrmax;
min_clr = (min_rel_mu_tmp-min_rel_mu)/range_mu*clrmax;
clr = [linspace(min_clr(1),max_clr(1),n_points);linspace(min_clr(2),max_clr(2),n_points);linspace(min_clr(3),max_clr(3),n_points)]';
colormap(ax2,clr);
colorbar;
set(ax2,'FontSize',12,'FontName','Helvetica');

% stage 3
ax3 = subplot(2,3,3);
data3 = rel_res1(:,:,idx_3);
x = st3.*ones(a,length(idx_3));
y = permute(data3(:,1,:),[1 3 2]);
z1 = permute(data3(:,2,:),[1 3 2]);
surf(x,y,z1,'EdgeColor',[1 1 1],'FaceColor','interp');
xlim([min(st3),max(st3)]);
ylim([1,1.5]);
zlim([1,1.7]);
set(gca,'xscale','log');
set(ax3,'xtick',[min(st3),max(st3)]);
set(ax3,'ztick',[1,1.5]);
view(az, el);
max_rel_mu_tmp = max(max(data3(:,2,:)));
min_rel_mu_tmp = min(min(data3(:,2,:)));
max_clr = (max_rel_mu_tmp-min_rel_mu)/range_mu*clrmax;
min_clr = (min_rel_mu_tmp-min_rel_mu)/range_mu*clrmax;
clr = [linspace(min_clr(1),max_clr(1),n_points);linspace(min_clr(2),max_clr(2),n_points);linspace(min_clr(3),max_clr(3),n_points)]';
colormap(ax3,clr);
colorbar;
set(ax3,'FontSize',12,'FontName','Helvetica');

% changing glucose transporter
% stage 1
ax4 = subplot(2,3,4);
data4 = rel_res2(:,:,idx_1);
x = st1.*ones(a,length(idx_1));
y = permute(data4(:,1,:),[1 3 2]);
z1 = permute(data4(:,2,:),[1 3 2]);
surf(x,y,z1,'EdgeColor',[1 1 1],'FaceColor','interp');
xlim([min(st1),max(st1)]);
ylim([1,1.5]);
zlim([1,1.7]);
set(ax4,'xtick',[min(st1),max(st1)]);
set(ax4,'ztick',[1,1.5]);
view(az, el);
max_rel_mu_tmp = max(max(data4(:,2,:)));
min_rel_mu_tmp = min(min(data4(:,2,:)));
max_clr = (max_rel_mu_tmp-min_rel_mu)/range_mu*clrmax;
min_clr = (min_rel_mu_tmp-min_rel_mu)/range_mu*clrmax;
clr = [linspace(min_clr(1),max_clr(1),n_points);linspace(min_clr(2),max_clr(2),n_points);linspace(min_clr(3),max_clr(3),n_points)]';
colormap(ax4,clr);
colorbar;
set(ax4,'FontSize',12,'FontName','Helvetica');

% stage 2
ax5 = subplot(2,3,5);
data5 = rel_res2(:,:,idx_2);
x = st2.*ones(a,length(idx_2));
y = permute(data5(:,1,:),[1 3 2]);
z1 = permute(data5(:,2,:),[1 3 2]);
surf(x,y,z1,'EdgeColor',[1 1 1],'FaceColor','interp');
xlim([min(st2),max(st2)]);
ylim([1,1.5]);
zlim([1,1.7]);
set(ax5,'xtick',[min(st2),max(st2)]);
set(ax5,'ztick',[1,1.5]);
view(az, el);
max_rel_mu_tmp = max(max(data5(:,2,:)));
min_rel_mu_tmp = min(min(data5(:,2,:)));
max_clr = (max_rel_mu_tmp-min_rel_mu)/range_mu*clrmax;
min_clr = (min_rel_mu_tmp-min_rel_mu)/range_mu*clrmax;
clr = [linspace(min_clr(1),max_clr(1),n_points);linspace(min_clr(2),max_clr(2),n_points);linspace(min_clr(3),max_clr(3),n_points)]';
colormap(ax5,clr);
colorbar;
set(ax5,'FontSize',12,'FontName','Helvetica');

% stage 3
ax6 = subplot(2,3,6);
data6 = rel_res2(:,:,idx_3);
x = st3.*ones(a,length(idx_3));
y = permute(data6(:,1,:),[1 3 2]);
z1 = permute(data6(:,2,:),[1 3 2]);
surf(x,y,z1,'EdgeColor',[1 1 1],'FaceColor','interp');
xlim([min(st3),max(st3)]);
ylim([1,1.5]);
zlim([1,1.7]);
set(gca,'xscale','log');
set(ax6,'xtick',[min(st3),max(st3)]);
set(ax6,'ztick',[1,1.5]);
view(az, el);
max_rel_mu_tmp = max(max(data6(:,2,:)));
min_rel_mu_tmp = min(min(data6(:,2,:)));
max_clr = (max_rel_mu_tmp-min_rel_mu)/range_mu*clrmax;
min_clr = (min_rel_mu_tmp-min_rel_mu)/range_mu*clrmax;
clr = [linspace(min_clr(1),max_clr(1),n_points);linspace(min_clr(2),max_clr(2),n_points);linspace(min_clr(3),max_clr(3),n_points)]';
colormap(ax6,clr);
colorbar;
set(ax6,'FontSize',12,'FontName','Helvetica');

%% Figures for analysing yield
