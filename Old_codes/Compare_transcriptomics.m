%% Compare predicted and measured protein levels

% The predicted data are from 'predicted_transcriptomics.mat'.

% The measured data are from 'Transcriptomics_Goel_2015.xlsx'.

%% Process predicted data
load('predicted_transcriptomics.mat');

% Without saturation factor
% Calculate mean and remove the data with coefficient of variation > 0.5
tmpdata = predicted.without_sf;
tmpdata(tmpdata == 0) = nan;
tmpdata(tmpdata < 1e-13) = nan;%remove extremely low values
ave015 = mean(tmpdata(:,1:3),2);
ave030 = mean(tmpdata(:,4:6),2);
ave050 = mean(tmpdata(:,8:9),2);
ave060 = tmpdata(:,10);
sd015 = std(tmpdata(:,1:3),0,2);
sd030 = std(tmpdata(:,4:6),0,2);
sd050 = std(tmpdata(:,8:9),0,2);
cv015 = sd015 ./ ave015;
cv030 = sd030 ./ ave030;
cv050 = sd050 ./ ave050;
ave015(cv015 > 0.5) = nan;
ave030(cv030 > 0.5) = nan;
ave050(cv050 > 0.5) = nan;
% Calculate log2foldchange using D=0.15 as reference
predicted.without_sf_rel = [log2(ave030./ave015)  log2(ave050./ave015) log2(ave060./ave015)];
clear ave015 ave030 ave050 ave060 cv015 cv030 cv050 tmpdata sd015 sd030 sd050;

% With saturation factor
% Calculate mean and remove the data with coefficient of variation > 0.5
tmpdata = predicted.with_sf;
tmpdata(tmpdata == 0) = nan;
tmpdata(tmpdata < 1e-13) = nan;%remove extremely low values
ave015 = mean(tmpdata(:,1:3),2);
ave030 = mean(tmpdata(:,4:6),2);
ave050 = mean(tmpdata(:,8:9),2);
ave060 = tmpdata(:,10);
sd015 = std(tmpdata(:,1:3),0,2);
sd030 = std(tmpdata(:,4:6),0,2);
sd050 = std(tmpdata(:,8:9),0,2);
cv015 = sd015 ./ ave015;
cv030 = sd030 ./ ave030;
cv050 = sd050 ./ ave050;
ave015(cv015 > 0.5) = nan;
ave030(cv030 > 0.5) = nan;
ave050(cv050 > 0.5) = nan;
% Calculate log2foldchange using D=0.15 as reference
predicted.with_sf_rel = [log2(ave030./ave015)  log2(ave050./ave015) log2(ave060./ave015)];
clear ave015 ave030 ave050 ave060 cv015 cv030 cv050 tmpdata sd015 sd030 sd050;

%% Extract measured data
[num, txt, ~] = xlsread('Transcriptomics_Goel_2015.xlsx','Processed_data_p0.05');

measured = struct();
measured.geneID = txt(2:end,1);
measured.rel = num(:,1:3);
clear num txt;

%% Compare processed datasets

combined = struct();
combined.geneID = measured.geneID(ismember(measured.geneID,predicted.geneID));
combined.rel = zeros(length(combined.geneID),9);
for i = 1:length(combined.geneID)
    geneid = combined.geneID(i);
    idx_measured = ismember(measured.geneID,geneid);
    measuredtmp = measured.rel(idx_measured,:);
    idx_predicted = ismember(predicted.geneID,geneid);
    predictedtmp_withoutsf = predicted.without_sf_rel(idx_predicted,:);
    predictedtmp_withsf = predicted.with_sf_rel(idx_predicted,:);
    combined.rel(i,:) = [measuredtmp predictedtmp_withoutsf predictedtmp_withsf];
end
clear i idx_measured idx_predicted measuredtmp predictedtmp_withoutsf predictedtmp_withsf protid;

measured_1 = combined.rel(:,1);
measured_2 = combined.rel(:,2);
measured_3 = combined.rel(:,3);
without_sf_1 = combined.rel(:,4);
without_sf_2 = combined.rel(:,5);
without_sf_3 = combined.rel(:,6);
with_sf_1 = combined.rel(:,7);
with_sf_2 = combined.rel(:,8);
with_sf_3 = combined.rel(:,9);
geneID_1 = cellfun(@(x) strcat(x,'_030'),combined.geneID,'UniformOutput',false);
geneID_2 = cellfun(@(x) strcat(x,'_050'),combined.geneID,'UniformOutput',false);
geneID_3 = cellfun(@(x) strcat(x,'_060'),combined.geneID,'UniformOutput',false);

measured_all = [measured_1;measured_2;measured_3];
without_sf_all = [without_sf_1;without_sf_2;without_sf_3];
with_sf_all = [with_sf_1;with_sf_2;with_sf_3];
geneID_all = [geneID_1;geneID_2;geneID_3];
clear measured_1 measured_2 measured_3;
clear without_sf_1 without_sf_2 without_sf_3;
clear with_sf_1 with_sf_2 with_sf_3;
clear geneID_1 geneID_2 geneID_3;

% Collect non-nan data and Calculate RMSE and R2
without_sf = struct();
I = ~isnan(measured_all) & ~isnan(without_sf_all);
without_sf.geneID = geneID_all(I);
without_sf.exp = measured_all(I);
without_sf.pred =  without_sf_all(I);
without_sf.rmse = sqrt(sum((without_sf.exp(:)-without_sf.pred(:)).^2)/numel(without_sf.exp));
[R,P] = corrcoef(without_sf.exp,without_sf.pred);
without_sf.p = P(1,2);
without_sf.R = round(R(1,2),3);
without_sf.R2 = round(R(1,2)^2,3);
without_sf.N = sum(I);
clear I R P;

with_sf = struct();
I = ~isnan(measured_all) & ~isnan(with_sf_all); 
with_sf.geneID = geneID_all(I);
with_sf.exp = measured_all(I);
with_sf.pred = with_sf_all(I);
with_sf.rmse = sqrt(sum((with_sf.exp(:)-with_sf.pred(:)).^2)/numel(with_sf.exp));
[R,P] = corrcoef(with_sf.exp,with_sf.pred);
with_sf.p = P(1,2);
with_sf.R = round(R(1,2),3);
with_sf.R2 = round(R(1,2)^2,3);
with_sf.N = sum(I);
clear I R P;

clear measured_all without_sf_all with_sf_all geneID_all;

%% Filter with metabolic proteins
% load('Gene_list.mat');
% I = contains(without_sf.geneID,Gene_list.M_gene);
% without_sf.geneID = without_sf.geneID(I);
% without_sf.exp = without_sf.exp(I);
% without_sf.pred =  without_sf.pred(I);
% without_sf.rmse = sqrt(sum((without_sf.exp(:)-without_sf.pred(:)).^2)/numel(without_sf.exp));
% [R,P] = corrcoef(without_sf.exp,without_sf.pred);
% without_sf.p = P(1,2);
% without_sf.R = round(R(1,2),3);
% without_sf.R2 = round(R(1,2)^2,3);
% without_sf.N = sum(I);
% clear I R P;
% 
% I = contains(with_sf.geneID,Gene_list.M_gene);
% with_sf.geneID = with_sf.geneID(I);
% with_sf.exp = with_sf.exp(I);
% with_sf.pred = with_sf.pred(I);
% with_sf.rmse = sqrt(sum((with_sf.exp(:)-with_sf.pred(:)).^2)/numel(with_sf.exp));
% [R,P] = corrcoef(with_sf.exp,with_sf.pred);
% with_sf.p = P(1,2);
% with_sf.R = round(R(1,2),3);
% with_sf.R2 = round(R(1,2)^2,3);
% with_sf.N = sum(I);
% clear I R P;

%% 
save('T_without_sf.mat','without_sf');
save('T_with_sf.mat','with_sf');

%% Figures for overview analyses
figure();
hold on;
box on;
x = without_sf.exp;
y = without_sf.pred;
p = polyfit(x,y,1);
f = polyval(p,x);
h = plot([min(x),max(x)],[min(f),max(f)],'-','Color',[69,117,180]/256,'LineWidth',5);
h.Color(4) = 0.5;
clear x y p f h;
scatter(without_sf.exp,without_sf.pred,8,'MarkerFaceColor',[69,117,180]/256,'MarkerEdgeColor',[69,117,180]/256,'MarkerFaceAlpha',1,'MarkerEdgeAlpha',0);
xlim([-14.9 14.9]);
ylim([-14.9 14.9]);
line([-14.9 14.9],[-14.9 14.9],'Color','k','LineStyle',':');
set(gca,'FontSize',12,'FontName','Helvetica');
title('without sf','FontSize',14,'FontName','Helvetica','fontweight','bold');
xlabel('Measured log_2FC','FontSize',14,'FontName','Helvetica');
ylabel('Predicted log_2FC','FontSize',14,'FontName','Helvetica');
str_rmse = ['RMSE = ',num2str(round(without_sf.rmse,2))];
str_r = ['R = ',num2str(round(without_sf.R,2))];
% str_p = ['p = ',num2str(without_sf.p)];
% str_n = ['N = ',num2str(without_sf.N)];
text(-1,-4.5,str_rmse,'FontSize',13,'FontName','Helvetica','Color',[69,117,180]/256);
text(-5,4.5,str_r,'FontSize',13,'FontName','Helvetica','Color',[69,117,180]/256);
% text(-0.5,-3,str_p,'FontSize',13,'FontName','Helvetica','Color',[69,117,180]/256);
% text(-0.5,0,str_n,'FontSize',13,'FontName','Helvetica','Color',[69,117,180]/256);
clear str_rmse str_r str_p str_n;
set(gcf,'position',[0 0 180 180]);
set(gca,'position',[0.23 0.23 0.75 0.75]);

figure();
hold on;
box on;
x = with_sf.exp;
y = with_sf.pred;
p = polyfit(x,y,1);
f = polyval(p,x);
h = plot([min(x),max(x)],[min(f),max(f)],'-','Color',[215,48,39]/256,'LineWidth',5);
h.Color(4) = 0.5;
clear x y p f h;
scatter(with_sf.exp,with_sf.pred,8,'MarkerFaceColor',[215,48,39]/256,'MarkerEdgeColor',[215,48,39]/256,'MarkerFaceAlpha',1,'MarkerEdgeAlpha',0);
xlim([-14.9 14.9]);
ylim([-14.9 14.9]);
line([-14.9 14.9],[-14.9 14.9],'Color','k','LineStyle',':');
set(gca,'FontSize',12,'FontName','Helvetica');
title('with sf','FontSize',14,'FontName','Helvetica','fontweight','bold');
xlabel('Measured log_2FC','FontSize',14,'FontName','Helvetica');
ylabel('Predicted log_2FC','FontSize',14,'FontName','Helvetica');
str_rmse = ['RMSE = ',num2str(round(with_sf.rmse,2))];
str_r = ['R = ',num2str(round(with_sf.R,2))];
% str_p = ['p = ',num2str(with_sf.p)];
% str_n = ['N = ',num2str(with_sf.N)];
text(-1,-4.5,str_rmse,'FontSize',13,'FontName','Helvetica','Color',[215,48,39]/256);
text(-5,4.5,str_r,'FontSize',13,'FontName','Helvetica','Color',[215,48,39]/256);
% text(-0.5,-3,str_p,'FontSize',13,'FontName','Helvetica','Color',[215,48,39]/256);
% text(-0.5,0,str_n,'FontSize',13,'FontName','Helvetica','Color',[215,48,39]/256);
clear str_rmse str_r str_p str_n;
set(gcf,'position',[250 0 180 180]);
set(gca,'position',[0.23 0.23 0.75 0.75]);

