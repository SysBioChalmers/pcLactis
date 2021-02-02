load('param_list.mat');
param_list = cellfun(@(x) strrep(x,'_Enzyme_c',''),param_list,'UniformOutput',false);
param_list = cellfun(@(x) strrep(x,'_','-'),param_list,'UniformOutput',false);
load('other_param_list.mat');
other_param_list = cellfun(@(x) strrep(x,'_',' '),other_param_list,'UniformOutput',false);

tot_param_list = [param_list;other_param_list];

load('mu_ref_1000.mat');
load('mu_ref_1747.mat');

load('mu_list_kcat_sens_dn_1000.mat');
load('mu_list_others_sens_dn_1000.mat');
mu_list_dn_1000 = [mu_list_kcat_sens_dn_1000;mu_list_others_sens_dn_1000];

load('mu_list_kcat_sens_up_1000.mat');
load('mu_list_others_sens_up_1000.mat');
mu_list_up_1000 = [mu_list_kcat_sens_up_1000;mu_list_others_sens_up_1000];

load('mu_list_kcat_sens_up_1747.mat');
load('mu_list_others_sens_up_1747.mat');
mu_list_up_1747 = [mu_list_kcat_sens_up_1747;mu_list_others_sens_up_1747];

load('mu_list_kcat_sens_dn_1747.mat');
load('mu_list_others_sens_dn_1747.mat');
mu_list_dn_1747 = [mu_list_kcat_sens_dn_1747;mu_list_others_sens_dn_1747];

rel_change_dn_1000 = mu_list_dn_1000/mu_ref_1000-1;
rel_change_up_1000 = mu_list_up_1000/mu_ref_1000-1;
rel_change_dn_1747 = mu_list_dn_1747/mu_ref_1747-1;
rel_change_up_1747 = mu_list_up_1747/mu_ref_1747-1;

figure();
subplot(1,2,1); 
data_raw = rel_change_dn_1747;
label_raw = tot_param_list;
[~, idx] = sort(abs(data_raw),'ascend');
data = data_raw(idx)*100;
label = label_raw(idx);
select_idx = 816:845;
h = barh(1:length(select_idx),data(select_idx),0.4,'FaceColor',[33,102,172]/255,'EdgeColor',[33,102,172]/255,'LineWidth',0.5);
xlim([-100 100]);
set(gca,'YTick',1:1:length(select_idx));
set(gca,'YTickLabel',label(select_idx));
set(gca,'FontSize',7,'FontName','Helvetica');
xlabel('Relative change in growth (%)','FontSize',7,'FontName','Helvetica','Color','k');
title('2-fold decrease','FontSize',7,'FontName','Helvetica','Color','k');
subplot(1,2,2); 
data_raw = rel_change_up_1747;
label_raw = tot_param_list;
[~, idx] = sort(abs(data_raw),'ascend');
data = data_raw(idx)*100;
label = label_raw(idx);
select_idx = 816:845;
h = barh(1:length(select_idx),data(select_idx),0.4,'FaceColor',[178,24,43]/255,'EdgeColor',[178,24,43]/255,'LineWidth',0.5);
xlim([-100 100]);
set(gca,'YTick',1:1:length(select_idx));
set(gca,'YTickLabel',label(select_idx));
set(gca,'FontSize',7,'FontName','Helvetica');
xlabel('Relative change in growth (%)','FontSize',7,'FontName','Helvetica','Color','k');
title('2-fold increase','FontSize',7,'FontName','Helvetica','Color','k');
set(gcf,'position',[450 50 550 380]);


figure();
subplot(1,2,1); 
data_raw = rel_change_dn_1000;
label_raw = tot_param_list;
[~, idx] = sort(abs(data_raw),'ascend');
data = data_raw(idx)*100;
label = label_raw(idx);
select_idx = 816:845;
h = barh(1:length(select_idx),data(select_idx),0.4,'FaceColor',[33,102,172]/255,'EdgeColor',[33,102,172]/255,'LineWidth',0.5);
xlim([-100 100]);
set(gca,'YTick',1:1:length(select_idx));
set(gca,'YTickLabel',label(select_idx));
set(gca,'FontSize',7,'FontName','Helvetica');
xlabel('Relative change in growth (%)','FontSize',7,'FontName','Helvetica','Color','k');
title('2-fold decrease','FontSize',7,'FontName','Helvetica','Color','k');
subplot(1,2,2); 
data_raw = rel_change_up_1000;
label_raw = tot_param_list;
[~, idx] = sort(abs(data_raw),'ascend');
data = data_raw(idx)*100;
label = label_raw(idx);
select_idx = 816:845;
h = barh(1:length(select_idx),data(select_idx),0.4,'FaceColor',[178,24,43]/255,'EdgeColor',[178,24,43]/255,'LineWidth',0.5);
xlim([-100 100]);
set(gca,'YTick',1:1:length(select_idx));
set(gca,'YTickLabel',label(select_idx));
set(gca,'FontSize',7,'FontName','Helvetica');
xlabel('Relative change in growth (%)','FontSize',7,'FontName','Helvetica','Color','k');
title('2-fold increase','FontSize',7,'FontName','Helvetica','Color','k');
set(gcf,'position',[450 550 550 380]);








