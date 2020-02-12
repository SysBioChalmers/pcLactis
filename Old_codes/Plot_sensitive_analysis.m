%% Plot sensitivity of glucose transporter and modeled proteome.

% Data import
load('Samp_result.mat');
load('Sagt_result.mat');

res1 = res_mp;
res2 = res_gt;

glc_list = [2 4 6 8 10 20 40 60 80 100 200 1000 10000 100000 1000000];%unit: /uM
idx_1 = 1:5;
idx_2 = 6:10;
idx_3 = 11:length(glc_list);
st1 = glc_list(idx_1);%unit: /uM
st2 = glc_list(idx_2);%unit: /uM
st3 = glc_list(idx_3);%unit: /uM

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

% Figures
figure('position',[0 0 1200 600]);
az = -20;
el = 20;
clrmax = [227,26,28]/255;
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
ylim([1,1.05]);
zlim([1,1.1]);
set(ax1,'xtick',[min(st1),max(st1)]);
set(ax1,'ztick',[1,1.1]);
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
ylim([1,1.05]);
zlim([1,1.1]);
set(ax2,'xtick',[min(st2),max(st2)]);
set(ax2,'ztick',[1,1.1]);
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
ylim([1,1.05]);
zlim([1,1.1]);
set(gca,'xscale','log');
set(ax3,'xtick',[min(st3),max(st3)]);
set(ax3,'ztick',[1,1.1]);
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
ylim([1,1.05]);
zlim([1,1.1]);
set(ax4,'xtick',[min(st1),max(st1)]);
set(ax4,'ztick',[1,1.1]);
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
ylim([1,1.05]);
zlim([1,1.1]);
set(ax5,'xtick',[min(st2),max(st2)]);
set(ax5,'ztick',[1,1.1]);
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
ylim([1,1.05]);
zlim([1,1.1]);
set(gca,'xscale','log');
set(ax6,'xtick',[min(st3),max(st3)]);
set(ax6,'ztick',[1,1.1]);
view(az, el);
max_rel_mu_tmp = max(max(data6(:,2,:)));
min_rel_mu_tmp = min(min(data6(:,2,:)));
max_clr = (max_rel_mu_tmp-min_rel_mu)/range_mu*clrmax;
min_clr = (min_rel_mu_tmp-min_rel_mu)/range_mu*clrmax;
clr = [linspace(min_clr(1),max_clr(1),n_points);linspace(min_clr(2),max_clr(2),n_points);linspace(min_clr(3),max_clr(3),n_points)]';
colormap(ax6,clr);
colorbar;
set(ax6,'FontSize',12,'FontName','Helvetica');


