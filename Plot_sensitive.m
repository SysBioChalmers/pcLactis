% Data import
load('Samp_result.mat');
load('Sagt_result.mat');

constraint_list_mp = res_mp(:,1,1);
constraint_list_mp = constraint_list_mp/constraint_list_mp(1);
mu_list_mp = permute(res_mp(:,2,:),[1 3 2]);
mu_list_mp = mu_list_mp./mu_list_mp(1,:);

constraint_list_gt = res_gt(:,1,1);
constraint_list_gt = constraint_list_gt/constraint_list_gt(1);
mu_list_gt = permute(res_gt(:,2,:),[1 3 2]);
mu_list_gt = mu_list_gt./mu_list_gt(1,:);

glc_list = [2 5 10 20 50 100 1000 10000 100000];
glc_label = {'2' '5' '10' '20' '50' '10^2' '10^3' '10^4' '10^5'};

az = -22;
el = 38;
clrmap = jet;
[n_points, ~] = size(clrmap);

max_rel_mu_mp = max(max(mu_list_mp));
min_rel_mu_mp = min(min(mu_list_mp));

max_rel_mu_gt = max(max(mu_list_gt));
min_rel_mu_gt = min(min(mu_list_gt));

max_rel_mu = max([max_rel_mu_mp,max_rel_mu_gt]);
min_rel_mu = min([min_rel_mu_mp,min_rel_mu_gt]);
range_mu = max_rel_mu-0;

interval = range_mu/n_points;


ax1 = subplot(1,2,1);
b = bar3(constraint_list_mp,mu_list_mp,0.3);
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
ylim([0.95,1.55]);
zlim([0.95,1.45]);
set(ax1,'YDir','normal');
set(ax1,'XTickLabel',glc_label);
set(ax1,'ytick',1:0.1:1.5);
set(ax1,'FontSize',10,'FontName','Helvetica');
zlabel('Change in growth rate','FontSize',12,'FontName','Helvetica');
ylabel('Change in constraint','FontSize',12,'FontName','Helvetica');
xlabel('Extracellular glucose concentration (uM)','FontSize',12,'FontName','Helvetica');
title('Modeled proteome','FontSize',14,'FontName','Helvetica');
view(az, el);
firstpoint = 1 + round((min_rel_mu_mp-min_rel_mu)/interval);
lastpoint = n_points + round((max_rel_mu_mp-max_rel_mu)/interval);
clrmap_mp = clrmap(firstpoint:lastpoint,:);
colormap(ax1,clrmap_mp);

ax2 = subplot(1,2,2);
b = bar3(constraint_list_gt,mu_list_gt,0.3);
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
ylim([0.95,1.55]);
zlim([0.95,1.75]);
set(ax2,'YDir','normal');
set(ax2,'XTickLabel',glc_label);
set(ax2,'ytick',1:0.1:1.5);
set(ax2,'FontSize',10,'FontName','Helvetica');
zlabel('Change in growth rate','FontSize',12,'FontName','Helvetica');
ylabel('Change in constraint','FontSize',12,'FontName','Helvetica');
xlabel('Extracellular glucose concentration (uM)','FontSize',12,'FontName','Helvetica');
title('Glucose transporter','FontSize',14,'FontName','Helvetica');
view(az, el);

firstpoint = 1 + round((min_rel_mu_gt-min_rel_mu)/interval);
lastpoint = n_points + round((max_rel_mu_gt-max_rel_mu)/interval);
clrmap_gt = clrmap(firstpoint:lastpoint,:);
colormap(ax2,clrmap_gt);

