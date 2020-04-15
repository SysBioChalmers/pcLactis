%% Data for modeled proteome and glucose transporter

f_unmodeled_ref = 0.4;
f_transporter_ref = 0.0083;
factor = 1.01;
modeled_protein_list = (1-f_unmodeled_ref)*[1,factor];
glucose_transporter_list = f_transporter_ref*[1,factor];
clear factor;

load('Sglc_result.mat');
selected_points = [1:2:23,25,26];
mu_list = glc_conc_without_sf(selected_points,1);
mu_tmp = mat2str(mu_list);
mu_tmp = mu_tmp(2:end-1);
mu_label = strsplit(mu_tmp,';')';

glc_list = glc_conc_without_sf(selected_points,2);
clear mu_tmp glc_conc_without_sf selected_points;

load('Sagt_fluxes.mat');
load('pcLactis_Model.mat');
model = pcLactis_Model;
load('Info_enzyme.mat');
glucTransporter = CalculateGlucoseTransporter(model,fluxes_gt,Info_enzyme);
glcT = glucTransporter/0.46;
clear glucTransporter model pcLactis_Model Info_enzyme;
clear fluxes_gt;

%% Scaled data
load('Samp_result.mat');
load('Sagt_result.mat');

% (Δmu/Δconstraint)*(ref_constraint/ref_mu)
% modeled proteome
delta_mu = res_mp(2,:) - res_mp(1,:);
delta_constraint = modeled_protein_list(2) - modeled_protein_list(1);
ref_constraint = modeled_protein_list(1);
ref_mu = res_mp(1,:);
data_mp = (delta_mu/delta_constraint).*(ref_constraint./ref_mu);

% glucose transporter
delta_mu = res_gt(2,:) - res_gt(1,:);
delta_constraint = glcT(2:2:end) - glcT(1:2:end);
ref_constraint = glcT(1:2:end);
ref_mu = res_gt(1,:);
data_gt = (delta_mu./delta_constraint).*(ref_constraint./ref_mu);

%% Plot

figure('Name','1');
minclr = [254,235,226]/255;
maxclr = [197,27,138]/255;
tmp1 = linspace(minclr(1),maxclr(1),129)';
tmp2 = linspace(minclr(2),maxclr(2),129)';
tmp3 = linspace(minclr(3),maxclr(3),129)';
clrmap = [tmp1 tmp2 tmp3];
h = heatmap(mu_label,{'Glucose transporter';'Modelled proteome'},[data_gt;data_mp],'Colormap',clrmap,...
    'ColorMethod','count','CellLabelColor','none');
h.XLabel = 'Growth rate (/h)';
title('Sensitivity analysis');
set(h,'FontSize',6,'FontName','Helvetica');
h.FontColor = 'k';
set(gcf,'position',[100 800 300 100]);
set(gca,'position',[0.2 0.4 0.4 0.16]);

% figure('Name','2');
% minclr = [199,234,229]/255;
% maxclr = [1,102,94]/255;
% tmp1 = linspace(minclr(1),maxclr(1),129)';
% tmp2 = linspace(minclr(2),maxclr(2),129)';
% tmp3 = linspace(minclr(3),maxclr(3),129)';
% clrmap = [tmp1 tmp2 tmp3];
% h = heatmap({'z'},mu_label,data_mp','Colormap',clrmap,...
%     'ColorMethod','count','CellLabelColor','none');
% title('Modelled proteome');
% set(h,'FontSize',6,'FontName','Helvetica');
% set(gcf,'position',[200 600 100 200]);
% set(gca,'position',[0.4 0.1 0.12 0.8]);
% 
% figure('Name','3');
% minclr = [255,255,255]/255;
% maxclr = [8,81,156]/255;
% tmp1 = linspace(minclr(1),maxclr(1),129)';
% tmp2 = linspace(minclr(2),maxclr(2),129)';
% tmp3 = linspace(minclr(3),maxclr(3),129)';
% clrmap = [tmp1 tmp2 tmp3];
% h = heatmap({'z'},mu_label,data_gt','Colormap',clrmap,...
%     'ColorMethod','count','CellLabelColor','none');
% title('Glucose transporter');
% set(h,'FontSize',6,'FontName','Helvetica');
% set(gcf,'position',[500 600 100 200]);
% set(gca,'position',[0.4 0.1 0.12 0.8]);