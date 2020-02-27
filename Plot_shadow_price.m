%% Data import and process
load('Samp_result.mat');
load('Sagt_result.mat');

constraint_list_mp = res_mp(:,1,1);
mu_list_mp = permute(res_mp(:,2,:),[1 3 2]);

constraint_list_gt = res_gt(:,1,1);
mu_list_gt = permute(res_gt(:,2,:),[1 3 2]);

glc_list = [2 5 10 20 50 100 1000 10000 100000];
glc_label = {'2' '5' '10' '20' '50' '10^2' '10^3' '10^4' '10^5'};

load('Sagt_fluxes.mat');
load('pcLactis_Model.mat');
model = pcLactis_Model;
load('Info_enzyme.mat');
glucTransporter = CalculateGlucoseTransporter(model,fluxes_gt,Info_enzyme);
[a, b] = size(mu_list_gt);
glcT = reshape(glucTransporter,a,b)/0.46;
clear glucTransporter model pcLactis_Model Info_enzyme;
clear res_gt res_mp fluxes_gt;

%% Scaled data
% (ref_constraint/ref_mu)*(Δmu/Δconstraint)
ref_mu_gt = mu_list_gt(1,:);
ref_gt_gt = glcT(1,:);
ratio_ref_gt = ref_gt_gt./ref_mu_gt;

ref_mu_mp = mu_list_mp(1,:);
ratio_ref_mp = constraint_list_mp(1)./ref_mu_mp;

data_gt = zeros(a-1, b);
for i = 2:a
    delta_mu = round((mu_list_gt(i,:) - ref_mu_gt),6);
    delta_constraint = round((glcT(i,:) - ref_gt_gt),6);
    ratio_delta = delta_mu./delta_constraint;
    scaled_value = ratio_ref_gt.*ratio_delta;
    data_gt(i-1,:) = scaled_value;
end

data_mp = zeros(a-1, b);
for i = 2:a
    delta_mu = round((mu_list_mp(i,:) - ref_mu_mp),6);
    delta_constraint = round((constraint_list_mp(i) - constraint_list_mp(1)),6);
    ratio_delta = delta_mu/delta_constraint;
    scaled_value = ratio_ref_mp.*ratio_delta;
    data_mp(i-1,:) = scaled_value;
end
clear a b;
clear delta_mu delta_constraint ratio_delta scaled_value i;

%% Plot

figure('Name','1');
minclr = [222,235,247]/255;
maxclr = [8,81,156]/255;
tmp1 = linspace(minclr(1),maxclr(1),129)';
tmp2 = linspace(minclr(2),maxclr(2),129)';
tmp3 = linspace(minclr(3),maxclr(3),129)';
clrmap = [tmp1 tmp2 tmp3];
h = heatmap(glc_label,{'1.1x','1.2x','1.3x','1.4x','1.5x'},data_gt,'Colormap',clrmap,'ColorMethod','count','CellLabelColor','none');
set(h,'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 300 350 120]);
set(gca,'position',[0.1 0.2 0.68 0.78]);

figure('Name','2');
minclr = [229,245,224]/255;
maxclr = [0,109,44]/255;
tmp1 = linspace(minclr(1),maxclr(1),129)';
tmp2 = linspace(minclr(2),maxclr(2),129)';
tmp3 = linspace(minclr(3),maxclr(3),129)';
clrmap = [tmp1 tmp2 tmp3];
h = heatmap(glc_label,{'1.1x','1.2x','1.3x','1.4x','1.5x'},data_mp,'Colormap',clrmap,'ColorMethod','count','CellLabelColor','none');
set(h,'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 600 350 120]);
set(gca,'position',[0.1 0.2 0.68 0.78]);


