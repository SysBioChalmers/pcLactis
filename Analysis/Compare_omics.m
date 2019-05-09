%% Compare proteomics and transcriptomics


%% Load data

load('P_without_sf.mat');
load('P_with_sf.mat');
P_without_sf = without_sf;
P_with_sf = with_sf;
clear without_sf with_sf;

load('T_without_sf.mat');
load('T_with_sf.mat');
T_without_sf = without_sf;
T_with_sf = with_sf;
clear without_sf with_sf;

%% Compare
% without sf
P_id_full = P_without_sf.proteinID;
T_id_full = T_without_sf.geneID;
P_exp_full = P_without_sf.exp;
T_exp_full = T_without_sf.exp;
P_pred_full = P_without_sf.pred;
T_pred_full = T_without_sf.pred;

result = struct();
result.id = P_id_full(ismember(P_id_full,T_id_full));
result.P_exp = zeros(length(result.id),1);
result.T_exp = zeros(length(result.id),1);
result.P_pred = zeros(length(result.id),1);
result.T_pred = zeros(length(result.id),1);
for i = 1:length(result.id)
    id_tmp = result.id(i);
    result.P_exp(i) = P_exp_full(ismember(P_id_full,id_tmp));
    result.T_exp(i) = T_exp_full(ismember(T_id_full,id_tmp));
    result.P_pred(i) = P_pred_full(ismember(P_id_full,id_tmp));
    result.T_pred(i) = T_pred_full(ismember(T_id_full,id_tmp));
end

figure();
scatter(result.P_exp,result.T_exp);
xlim([-5.5 5.5]);
ylim([-5.5 5.5]);
line([-5.5 5.5],[-5.5 5.5],'Color','k','LineStyle',':');
set(gca,'FontSize',12,'FontName','Helvetica');
title('without sf exp','FontSize',14,'FontName','Helvetica','fontweight','bold');
xlabel('log_2FC protein','FontSize',14,'FontName','Helvetica');
ylabel('log_2FC transcript','FontSize',14,'FontName','Helvetica');
set(gcf,'position',[0 0 300 300]);
set(gca,'position',[0.15 0.15 0.75 0.75]);
figure();
scatter(result.P_pred,result.T_pred);
xlim([-5.5 5.5]);
ylim([-5.5 5.5]);
line([-5.5 5.5],[-5.5 5.5],'Color','k','LineStyle',':');
set(gca,'FontSize',12,'FontName','Helvetica');
title('without sf pred','FontSize',14,'FontName','Helvetica','fontweight','bold');
xlabel('log_2FC protein','FontSize',14,'FontName','Helvetica');
ylabel('log_2FC transcript','FontSize',14,'FontName','Helvetica');
set(gcf,'position',[300 0 300 300]);
set(gca,'position',[0.15 0.15 0.75 0.75]);

% with sf

P_id_full = P_with_sf.proteinID;
T_id_full = T_with_sf.geneID;
P_exp_full = P_with_sf.exp;
T_exp_full = T_with_sf.exp;
P_pred_full = P_with_sf.pred;
T_pred_full = T_with_sf.pred;

result = struct();
result.id = P_id_full(ismember(P_id_full,T_id_full));
result.P_exp = zeros(length(result.id),1);
result.T_exp = zeros(length(result.id),1);
result.P_pred = zeros(length(result.id),1);
result.T_pred = zeros(length(result.id),1);
for i = 1:length(result.id)
    id_tmp = result.id(i);
    result.P_exp(i) = P_exp_full(ismember(P_id_full,id_tmp));
    result.T_exp(i) = T_exp_full(ismember(T_id_full,id_tmp));
    result.P_pred(i) = P_pred_full(ismember(P_id_full,id_tmp));
    result.T_pred(i) = T_pred_full(ismember(T_id_full,id_tmp));
end

figure();
scatter(result.P_exp,result.T_exp);
xlim([-5.5 5.5]);
ylim([-5.5 5.5]);
line([-5.5 5.5],[-5.5 5.5],'Color','k','LineStyle',':');
set(gca,'FontSize',12,'FontName','Helvetica');
title('with sf exp','FontSize',14,'FontName','Helvetica','fontweight','bold');
xlabel('log_2FC protein','FontSize',14,'FontName','Helvetica');
ylabel('log_2FC transcript','FontSize',14,'FontName','Helvetica');
set(gcf,'position',[600 0 300 300]);
set(gca,'position',[0.15 0.15 0.75 0.75]);

figure();
scatter(result.P_pred,result.T_pred);
xlim([-5.5 5.5]);
ylim([-5.5 5.5]);
line([-5.5 5.5],[-5.5 5.5],'Color','k','LineStyle',':');
set(gca,'FontSize',12,'FontName','Helvetica');
title('with sf pred','FontSize',14,'FontName','Helvetica','fontweight','bold');
xlabel('log_2FC protein','FontSize',14,'FontName','Helvetica');
ylabel('log_2FC transcript','FontSize',14,'FontName','Helvetica');
set(gcf,'position',[900 0 300 300]);
set(gca,'position',[0.15 0.15 0.75 0.75]);

clear;