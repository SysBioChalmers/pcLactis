%% Compare predicted and measured protein levels

% The predicted data are from 'predicted.mat'.

% The measured data are from 'Proteomics_Goel_2015.xlsx'.

%% Process predicted data
load('predicted_proteomics1.mat');

% Without saturation factor
% Calculate mean and remove the data with coefficient of variation > 0.5
tmpdata = predicted.without_sf;
tmpdata(tmpdata == 0) = nan;
tmpdata(tmpdata < 1e-10) = nan;%remove extremely low values
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
tmpdata(tmpdata < 1e-10) = nan;%remove extremely low values
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

%% Process measured data
[num, txt, ~] = xlsread('Proteomics_Goel_2015.xlsx','Processed_data_p0.05');
idlist = txt(2:end,1);
unqid = unique(idlist);
locatelist = txt(2:end,2);

measured = struct();
measured.proteinID = [];
measured.rel = [];

load('mem_proteins.mat');
for i = 1:length(unqid)
    protid = unqid(i);
    tmp_rel = num(ismember(idlist,protid),1:3);
    if ismember(protid,mem_proteins)
        locateinfo = 'mem';
    else
        locateinfo = 'sol';
    end
    if ismember(locateinfo,locatelist(ismember(idlist,protid)))
        idxtmp = ismember(locatelist(ismember(idlist,protid)),locateinfo);
        rel = tmp_rel(idxtmp,:);
        measured.proteinID = [measured.proteinID;protid];
        measured.rel = [measured.rel;rel];
    end
end
clear i idlist idxtmp locateinfo locatelist num protid rel tmp_rel txt unqid mem_proteins;

%% Compare processed datasets

combined = struct();
combined.proteinID = measured.proteinID(ismember(measured.proteinID,predicted.proteinID));
combined.rel = zeros(length(combined.proteinID),9);
for i = 1:length(combined.proteinID)
    protid = combined.proteinID(i);
    idx_measured = ismember(measured.proteinID,protid);
    measuredtmp = measured.rel(idx_measured,:);
    idx_predicted = ismember(predicted.proteinID,protid);
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
proteinid_1 = cellfun(@(x) strcat(x,'_030'),combined.proteinID,'UniformOutput',false);
proteinid_2 = cellfun(@(x) strcat(x,'_050'),combined.proteinID,'UniformOutput',false);
proteinid_3 = cellfun(@(x) strcat(x,'_060'),combined.proteinID,'UniformOutput',false);

measured_all = [measured_1;measured_2;measured_3];
without_sf_all = [without_sf_1;without_sf_2;without_sf_3];
with_sf_all = [with_sf_1;with_sf_2;with_sf_3];
proteinid_all = [proteinid_1;proteinid_2;proteinid_3];
clear measured_1 measured_2 measured_3;
clear without_sf_1 without_sf_2 without_sf_3;
clear with_sf_1 with_sf_2 with_sf_3;
clear proteinid_1 proteinid_2 proteinid_3;

% Collect non-nan data and Calculate RMSE and R2
without_sf = struct();
I = ~isnan(measured_all) & ~isnan(without_sf_all);
without_sf.proteinID = proteinid_all(I);
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
with_sf.proteinID = proteinid_all(I);
with_sf.exp = measured_all(I);
with_sf.pred = with_sf_all(I);
with_sf.rmse = sqrt(sum((with_sf.exp(:)-with_sf.pred(:)).^2)/numel(with_sf.exp));
[R,P] = corrcoef(with_sf.exp,with_sf.pred);
with_sf.p = P(1,2);
with_sf.R = round(R(1,2),3);
with_sf.R2 = round(R(1,2)^2,3);
with_sf.N = sum(I);
clear I R P;

clear measured_all without_sf_all with_sf_all proteinid_all;

%% Filter with metabolic proteins
% load('Gene_list.mat');
% I = contains(without_sf.proteinID,Gene_list.M_gene);
% without_sf.proteinID = without_sf.proteinID(I);
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
% I = contains(with_sf.proteinID,Gene_list.M_gene);
% with_sf.proteinID = with_sf.proteinID(I);
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
save('P_without_sf.mat','without_sf');
save('P_with_sf.mat','with_sf');

%% Figures for overview analyses
color1 = [77,77,77]/255;
color2 = [178,24,43]/255;
% color2 = [221,28,119]/255;

figure();
hold on;
box on;
x = without_sf.exp;
y = without_sf.pred;
p = polyfit(x,y,1);
f = polyval(p,x);
h = plot([min(x),max(x)],[min(f),max(f)],'-','Color',color1,'LineWidth',3);
h.Color(4) = 0.5;
clear x y p f h;
scatter(without_sf.exp,without_sf.pred,8,'MarkerFaceColor',color1,'MarkerEdgeColor',color1,'MarkerFaceAlpha',1,'MarkerEdgeAlpha',0);
xlim([-5.5 5.5]);
ylim([-5.5 5.5]);
line([-5.5 5.5],[-5.5 5.5],'Color','k','LineStyle',':');
set(gca,'FontSize',12,'FontName','Helvetica');
title('without sf','FontSize',14,'FontName','Helvetica','fontweight','bold');
xlabel('Measured log_2FC','FontSize',14,'FontName','Helvetica');
ylabel('Predicted log_2FC','FontSize',14,'FontName','Helvetica');
str_rmse = ['RMSE = ',num2str(round(without_sf.rmse,2))];
str_r = ['R = ',num2str(round(without_sf.R,2))];
% str_p = ['p = ',num2str(without_sf.p)];
% str_n = ['N = ',num2str(without_sf.N)];
text(-1,-4.5,str_rmse,'FontSize',13,'FontName','Helvetica','Color',color1);
text(-5,4.5,str_r,'FontSize',13,'FontName','Helvetica','Color',color1);
% text(-0.5,-3,str_p,'FontSize',13,'FontName','Helvetica','Color',color1);
% text(-0.5,0,str_n,'FontSize',13,'FontName','Helvetica','Color',color1);
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
h = plot([min(x),max(x)],[min(f),max(f)],'-','Color',color2,'LineWidth',3);
h.Color(4) = 0.5;
clear x y p f h;
scatter(with_sf.exp,with_sf.pred,8,'MarkerFaceColor',color2,'MarkerEdgeColor',color2,'MarkerFaceAlpha',1,'MarkerEdgeAlpha',0);
xlim([-5.5 5.5]);
ylim([-5.5 5.5]);
line([-5.5 5.5],[-5.5 5.5],'Color','k','LineStyle',':');
set(gca,'FontSize',12,'FontName','Helvetica');
title('with sf','FontSize',14,'FontName','Helvetica','fontweight','bold');
xlabel('Measured log_2FC','FontSize',14,'FontName','Helvetica');
ylabel('Predicted log_2FC','FontSize',14,'FontName','Helvetica');
str_rmse = ['RMSE = ',num2str(round(with_sf.rmse,2))];
str_r = ['R = ',num2str(round(with_sf.R,2))];
% str_p = ['p = ',num2str(with_sf.p)];
% str_n = ['N = ',num2str(with_sf.N)];
text(-1,-4.5,str_rmse,'FontSize',13,'FontName','Helvetica','Color',color2);
text(-5,4.5,str_r,'FontSize',13,'FontName','Helvetica','Color',color2);
% text(-0.5,-3,str_p,'FontSize',13,'FontName','Helvetica','Color',color2);
% text(-0.5,0,str_n,'FontSize',13,'FontName','Helvetica','Color',color2);
clear str_rmse str_r str_p str_n;
set(gcf,'position',[250 0 180 180]);
set(gca,'position',[0.23 0.23 0.75 0.75]);

%% Analyses for pathways
[~, txt1, ~] = xlsread('KEGG_pathway_llm.xlsx','GeneInPathway');
[~, txt2, ~] = xlsread('KEGG_pathway_llm.xlsx','PathwayID');
pathway_info = struct();
pathway_info.ID = txt2(:,1);
pathway_info.name = txt2(:,2);
pathway_info.relation_PW = txt1(:,1);
pathway_info.relation_prot = txt1(:,2);
clear txt1 txt2;

pathway_analysis = struct();
pathway_analysis.pathwayID = pathway_info.ID;
pathway_analysis.pathwayname = pathway_info.name;
n_tmp = length(pathway_analysis.pathwayname);
pathway_analysis.total_num_proteins = zeros(n_tmp,1);
pathway_analysis.withoutsf_proteins = zeros(n_tmp,1);
pathway_analysis.withsf_proteins = zeros(n_tmp,1);
pathway_analysis.withoutsf_fisherp = zeros(n_tmp,1);
pathway_analysis.withsf_fisherp = zeros(n_tmp,1);
pathway_analysis.withoutsf_rmse = zeros(n_tmp,1);
pathway_analysis.withsf_rmse = zeros(n_tmp,1);
pathway_analysis.withoutsf_R = zeros(n_tmp,1);
pathway_analysis.withsf_R = zeros(n_tmp,1);
pathway_analysis.withoutsf_p = zeros(n_tmp,1);
pathway_analysis.withsf_p = zeros(n_tmp,1);
clear n_tmp;

n_1 = length(without_sf.proteinID);
n_2 = length(with_sf.proteinID);
n_total = 2339; % number of protein genes in L. lactis, from NCBI database

for i = 1:length(pathway_analysis.pathwayID)
    pwid = pathway_analysis.pathwayID(i);
    idxpw = ismember(pathway_info.relation_PW,pwid);
    protlist = pathway_info.relation_prot(idxpw);
    num_tot_prot = length(protlist);
    
    idx1 = contains(without_sf.proteinID,protlist);
    proteins1_tmp = without_sf.proteinID(idx1);
    proteins1_tmp = cellfun(@(x) x(1:9),proteins1_tmp,'UniformOutput',false);
    proteins1_tmp = unique(proteins1_tmp);
    num_prot1 = length(proteins1_tmp);
    
    %Fisher Exact P-Value is used for enrichment analysis, p < 0.05
    [h1,fisherp1,~] = fishertest([num_prot1,num_tot_prot;n_1-num_prot1,n_total-num_tot_prot]);
    
	without_exp = without_sf.exp(idx1);
	without_pred = without_sf.pred(idx1);
    if h1
        rmse1 = sqrt(sum((without_exp(:)-without_pred(:)).^2)/numel(without_exp));
        [Rlist1,Plist1] = corrcoef(without_exp,without_pred);
        R1 = Rlist1(1,2);
        P1 = Plist1(1,2);
    else
        rmse1 = nan;
        R1 = nan;
        P1 = nan;
    end
    
    idx2 = contains(with_sf.proteinID,protlist);
    proteins2_tmp = with_sf.proteinID(idx2);
    proteins2_tmp = cellfun(@(x) x(1:9),proteins2_tmp,'UniformOutput',false);
    proteins2_tmp = unique(proteins2_tmp);
    num_prot2 = length(proteins2_tmp);
    
    %Fisher Exact P-Value is used for enrichment analysis, p < 0.05
    [h2,fisherp2,~] = fishertest([num_prot2,num_tot_prot;n_2-num_prot2,n_total-num_tot_prot]);

    with_exp = with_sf.exp(idx2);
    with_pred = with_sf.pred(idx2);
    if h2
        rmse2 = sqrt(sum((with_exp(:)-with_pred(:)).^2)/numel(with_exp));
        [Rlist2,Plist2] = corrcoef(with_exp,with_pred);
        R2 = Rlist2(1,2);
        P2 = Plist2(1,2);
    else
        rmse2 = nan;
        R2 = nan;
        P2 = nan;
    end

    pathway_analysis.total_num_proteins(i) = num_tot_prot;
    pathway_analysis.withoutsf_proteins(i) = num_prot1;
    pathway_analysis.withsf_proteins(i) = num_prot2;
    pathway_analysis.withoutsf_fisherp(i) = fisherp1;
    pathway_analysis.withsf_fisherp(i) = fisherp2;
    pathway_analysis.withoutsf_rmse(i) = rmse1;
    pathway_analysis.withsf_rmse(i) = rmse2;
    pathway_analysis.withoutsf_R(i) = R1;
    pathway_analysis.withsf_R(i) = R2;
    pathway_analysis.withoutsf_p(i) = P1;
    pathway_analysis.withsf_p(i) = P2;
end
clear i idx1 idx2 idxpw num_prot1 num_prot2 num_tot_prot P1 P2 R1 R2;
clear Plist1 Plist2 protlist pwid Rlist1 Rlist2 rmse1 rmse2;
clear with_exp with_pred without_exp without_pred proteins1_tmp proteins2_tmp;
clear fisherp1 fisherp2 h1 h2 n_1 n_2 n_total;
    
I = ~isnan(pathway_analysis.withoutsf_rmse) & ~isnan(pathway_analysis.withsf_rmse)...
     & ~isnan(pathway_analysis.withoutsf_R) & ~isnan(pathway_analysis.withsf_R);
common = struct();
common.pathway = pathway_analysis.pathwayname(I);
common.withoutsf_rmse = pathway_analysis.withoutsf_rmse(I);
common.withsf_rmse = pathway_analysis.withsf_rmse(I);
common.withoutsf_R = pathway_analysis.withoutsf_R(I);
common.withsf_R = pathway_analysis.withsf_R(I);
clear I;

figure();
cdata = [common.withoutsf_rmse common.withsf_rmse];
xvalues = {'Before','After'};
yvalues = common.pathway;
map = colormap(heatmap(xvalues,yvalues,cdata,'Colormap',bone));
map = sort(map,'descend');
h = heatmap(xvalues,yvalues,cdata,'Colormap',map,'CellLabelColor','none');
h.Title = 'RMSE';
set(gca,'FontSize',13,'FontName','Helvetica');
set(gcf,'position',[600 0 500 500]);
set(gca,'position',[0.7 0.1 0.1 0.6]);
clear cdata h map xvalues yvalues;

figure();
cdata = [common.withoutsf_R common.withsf_R];
xvalues = {'Before','After'};
yvalues = common.pathway;
map = colormap(heatmap(xvalues,yvalues,cdata,'Colormap',pink));
map = sort(map,'descend');
h = heatmap(xvalues,yvalues,cdata,'Colormap',map,'CellLabelColor','none');
h.Title = 'R';
set(gca,'FontSize',13,'FontName','Helvetica');
set(gcf,'position',[800 0 500 500]);
set(gca,'position',[0.7 0.1 0.1 0.6]);
clear cdata h map xvalues yvalues;
