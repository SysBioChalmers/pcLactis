%% Estimate saturation factors using experimental data.

% The codes below can generate the saturation factor mat file, named as
% sat_factor.mat

% Here should use M model since the pcLactis model has isozyme reactions.
load('M_Model.mat');
model = M_Model;
% But pcLactis model is also needed for getting info.
load('pcLactis_Model.mat');
model_pcLactis = pcLactis_Model;

clear M_Model pcLactis_Model;

%% Setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set objective funtion.
model = changeObjective(model,'R_M_ATPM',1);

% Set NGAM.
model=changeRxnBounds(model,'R_M_ATPM',0,'l');
model=changeRxnBounds(model,'R_M_ATPM',1000,'u');

% Block uptake direction of all the exchange reactions.
idx = cellfun(@(x) contains(x,'R_M_EX_'),model.rxns,'UniformOutput',false);
M_exchange_reactions = model.rxns(cell2mat(idx));
model = changeRxnBounds(model,M_exchange_reactions,0,'l');
clear idx M_exchange_reactions;

% Block some reactions in M part.
model = changeRxnBounds(model,'R_M_biomass_LLA_atpm',0,'b');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loop for calculating fluxes at dilution rate of 0.15 0.3 0.5 and 0.6.

[~, ~, exchange_raw] = xlsread('Exchange_reaction_setting.xlsx','Exp_bounds_GAM_NGAM');

Exchange_reactions = exchange_raw(2:end,1);
header = exchange_raw(1,1:21);

flux_distribution = zeros(length(model.rxns),0);

mu_list = [0.15 0.3 0.5 0.6];

for i = 1:length(mu_list)
    mu = mu_list(i);
    
    model_tmp = changeRxnBounds(model,'R_M_biomass_LLA',mu,'b');

    % Set bounds for some exchange reactions.
    idx = find(contains(header,num2str(mu)));
    replicate = length(idx)/2;
    
    for j = 1:replicate
        lb_idx = idx(2*j-1);
        ub_idx = idx(2*j);
        LB = cell2mat(exchange_raw(2:end,lb_idx));
        UB = cell2mat(exchange_raw(2:end,ub_idx));
        model_tmp = changeRxnBounds(model_tmp,Exchange_reactions,LB,'l');
        model_tmp = changeRxnBounds(model_tmp,Exchange_reactions,UB,'u');

        sol_tmp = optimizeCbModel(model_tmp,'max','one');
        
        [~, k] = size(flux_distribution);
        flux_distribution(:,k+1) = sol_tmp.x;
    end
end


%% Data processing (fluxes data)

flux_data = struct;
flux_data.value_tmp = [];
flux_data.rxnID = {};
flux_data.GPR = {};
flux_data.dir = [];

% Only keep rxns with non-zero fluxes at all the dilution rates. Besides, 
% remove rxns without GPRs.
for i = 1:length(flux_distribution)
    row = flux_distribution(i,:);
    if all(row) && ~isempty(model.grRules{i}) ...
       && ~strcmp(model.grRules{i},'dummy') && (all(row > 0) || all(row < 0))
        flux_data.value_tmp = [flux_data.value_tmp;row];
        flux_data.rxnID = [flux_data.rxnID;model.rxns(i)];
        flux_data.GPR = [flux_data.GPR;model.grRules(i)];
        if all(row > 0)
            flux_data.dir = [flux_data.dir;1];
        else
            flux_data.dir = [flux_data.dir;-1];
        end
    end
end

flux_data.mean = zeros(length(flux_data.rxnID),0);
flux_data.mean(:,1) = mean(flux_data.value_tmp(:,1:3),2);
flux_data.mean(:,2) = mean(flux_data.value_tmp(:,4:6),2);
flux_data.mean(:,3) = mean(flux_data.value_tmp(:,7:8),2);
flux_data.mean(:,4) = flux_data.value_tmp(:,9);

flux_data.logFC = zeros(length(flux_data.rxnID),0);
flux_data.logFC(:,1) = log2(flux_data.mean(:,2)./flux_data.mean(:,1));
flux_data.logFC(:,2) = log2(flux_data.mean(:,3)./flux_data.mean(:,1));
flux_data.logFC(:,3) = log2(flux_data.mean(:,4)./flux_data.mean(:,1));

%% Data processing (analysis with proteomics data)
[~, ~, proteomics_raw] = xlsread('Proteomics_Goel_2015','Processed_data');

cmp = struct;
cmp.protein_logFC = [];
cmp.RxnID = {};
cmp.flux_logFC = [];
cmp.ProteinID = {};
cmp.flux_dir = [];

for i = 1:length(flux_data.GPR)
    
    % determine if the rxn contains extracellular mets
    rxn_id = flux_data.rxnID{i};
    met_id = model.mets(model.S(:,strcmp(model.rxns,rxn_id)) ~= 0);
    num_met_e = sum(cell2mat(cellfun(@(x) strcmp(x(end),'e'),met_id,'UniformOutput',false)));
    if num_met_e == 0
        loc = 'sol';
    else
        loc = 'mem';
    end
    
    % compare flux and protein
    GPRs_tmp = flux_data.GPR{i};
    GPRs_tmp = strrep(GPRs_tmp,'(',' ');
    GPRs_tmp = strrep(GPRs_tmp,')',' ');
    GPRs_tmp = strrep(GPRs_tmp,'and',' ');
    GPRs_tmp = strrep(GPRs_tmp,'or',' ');
    GPRs_tmp = strtrim(GPRs_tmp);
    sub = unique(strsplit(GPRs_tmp));
    
    for j = 1:length(sub)
        sub_id = sub{j};
        
        if ~all(strcmp(proteomics_raw(:,1),sub_id) == 0)
            data_tmp = proteomics_raw(strcmp(proteomics_raw(:,1),sub_id),:);
            
            if sum(strcmp(data_tmp(:,3),loc)) > 0
                data_tmp_tmp = data_tmp(strcmp(data_tmp(:,3),loc),:);
                protein_value = cell2mat(data_tmp_tmp(1,4:6));
                
                cmp.protein_logFC = [cmp.protein_logFC;protein_value];
                cmp.RxnID = [cmp.RxnID;rxn_id];
                cmp.flux_logFC = [cmp.flux_logFC;flux_data.logFC(i,:)];
                cmp.ProteinID = [cmp.ProteinID;sub_id];
                cmp.flux_dir = [cmp.flux_dir;flux_data.dir(i)];
            end
        end
    end
end
                
% link rxnID in M model to pcLactis model
load('Info_enzyme.mat');

cmp_pcLactis = struct;
cmp_pcLactis.protein_logFC = [];
cmp_pcLactis.EnzymeID_pc = {};
cmp_pcLactis.flux_logFC = [];
cmp_pcLactis.ProteinID = {};

for i = 1:length(cmp.RxnID)
    rxn_id = cmp.RxnID{i};
    protein_id = cmp.ProteinID{i};

    gpr_list_pc = model_pcLactis.grRules(contains(model_pcLactis.rxns,rxn_id));
                
    if model.rev(strcmp(model.rxns,rxn_id)) == 1
        rxn_dir = cmp.flux_dir(i);
        if rxn_dir == 1
            dir_str = '_fwd_';
        elseif rxn_dir == -1
            dir_str = '_rvs_';
        end
        gpr_tmp = gpr_list_pc(contains(gpr_list_pc,dir_str));
    else
        gpr_tmp = gpr_list_pc;
    end
        
	for j = 1:length(gpr_tmp)
        gpr_id = gpr_tmp{j};
        
        if ~isempty(gpr_id)
            sub_list = Info_enzyme.subunit{strcmp(Info_enzyme.ID,gpr_id)};
            
            if sum(contains(sub_list,protein_id)) > 0
                cmp_pcLactis.EnzymeID_pc = [cmp_pcLactis.EnzymeID_pc;gpr_id];
                cmp_pcLactis.protein_logFC = [cmp_pcLactis.protein_logFC;cmp.protein_logFC(i,:)];
                cmp_pcLactis.flux_logFC = [cmp_pcLactis.flux_logFC;cmp.flux_logFC(i,:)];
                cmp_pcLactis.ProteinID = [cmp_pcLactis.ProteinID;protein_id];
            end
        end
	end
end

%% Estimate saturation factor.

% calculate correlation coefficient
cmp_pcLactis.corr_coef = [];
cmp_pcLactis.factor_ratio = [];

for i = 1:length(cmp_pcLactis.EnzymeID_pc)
    Ratio030015_V = 2^cmp_pcLactis.flux_logFC(i,1);
    Ratio050015_V = 2^cmp_pcLactis.flux_logFC(i,2);
    Ratio060015_V = 2^cmp_pcLactis.flux_logFC(i,3);

    Ratio030015_E = 2^cmp_pcLactis.protein_logFC(i,1);
    Ratio050015_E = 2^cmp_pcLactis.protein_logFC(i,2);
    Ratio060015_E = 2^cmp_pcLactis.protein_logFC(i,3);
    
    % V = factor * kcat * E, in which kcat is constant for each enzyme
    Ratio030015_Factor = Ratio030015_V / Ratio030015_E;
    Ratio050015_Factor = Ratio050015_V / Ratio050015_E;
    Ratio060015_Factor = Ratio060015_V / Ratio060015_E;
    
    ratio_factor = [Ratio030015_Factor Ratio050015_Factor Ratio060015_Factor];
    cmp_pcLactis.factor_ratio = [cmp_pcLactis.factor_ratio;ratio_factor];
    
    R_tmp = corrcoef([0.15 0.3 0.5 0.6],[1 ratio_factor]);
    R = R_tmp(1,2);
    cmp_pcLactis.corr_coef = [cmp_pcLactis.corr_coef;R];
end

% remove replicate with low abs corr coef

res_sat_factor = struct;
res_sat_factor.EnzymeID_pc = {};
res_sat_factor.ProteinID = {};
res_sat_factor.factor_ratio = [];
res_sat_factor.corr_coef = [];

tmp = unique(cmp_pcLactis.EnzymeID_pc);

for i = 1:length(tmp)
    enzyme_id = tmp{i};
    protein_tmp = cmp_pcLactis.ProteinID(strcmp(cmp_pcLactis.EnzymeID_pc,enzyme_id));
    factoratio_tmp = cmp_pcLactis.factor_ratio(strcmp(cmp_pcLactis.EnzymeID_pc,enzyme_id),:);
    coeff_tmp = cmp_pcLactis.corr_coef(strcmp(cmp_pcLactis.EnzymeID_pc,enzyme_id));
    [~, idx_max] = max(abs(coeff_tmp));
    
    res_sat_factor.EnzymeID_pc = [res_sat_factor.EnzymeID_pc;enzyme_id];
    res_sat_factor.ProteinID = [res_sat_factor.ProteinID;protein_tmp(idx_max)];
    res_sat_factor.factor_ratio = [res_sat_factor.factor_ratio;factoratio_tmp(idx_max,:)];
    res_sat_factor.corr_coef = [res_sat_factor.corr_coef;coeff_tmp(idx_max)];
end


%% Result processing.
% Remove enzymes with abs corr coef < 0.8
sat_factor = struct;
sat_factor.EnzymeID_pc = {};
sat_factor.ProteinID = {};
sat_factor.factor_ratio = [];
sat_factor.corr_coef = [];
sat_factor.factor = [];
sat_factor.poor_Enzyme = {};

idx = res_sat_factor.corr_coef >= 0.8;
sat_factor.EnzymeID_pc = res_sat_factor.EnzymeID_pc(idx);
sat_factor.ProteinID = res_sat_factor.ProteinID(idx);
sat_factor.factor_ratio = res_sat_factor.factor_ratio(idx,:);
sat_factor.corr_coef = res_sat_factor.corr_coef(idx);

% Enzymes with poor positive correlation (corr < 0.8)
idx_poor = res_sat_factor.corr_coef < 0.8;
sat_factor.poor_Enzyme = res_sat_factor.EnzymeID_pc(idx_poor);

% Determine linear equation with intercept of 0 for each enzyme, i.e., y=kx
if min(sat_factor.corr_coef) > 0
    sat_factor.factor(:,4) = ones(length(sat_factor.corr_coef),1);
    sat_factor.factor(:,1) = sat_factor.factor(:,4) ./ sat_factor.factor_ratio(:,3);
    sat_factor.factor(:,2) = sat_factor.factor(:,1) .* sat_factor.factor_ratio(:,1);
    sat_factor.factor(:,3) = sat_factor.factor(:,1) .* sat_factor.factor_ratio(:,2);
else
    error('should check the factors with negative correlation with mu');
end

sat_factor.k = [];
mu = [0.15;0.3;0.5;0.6];
for i = 1:length(sat_factor.EnzymeID_pc)
    y = transpose(sat_factor.factor(i,:));
    coef_tmp = mu\y;
    sat_factor.k = [sat_factor.k;coef_tmp];
end

cd Simulations/;
save('sat_factor.mat','sat_factor');
cd ../;

%% Figure distribution of all correlationss
figure();
hold on;
box on;
[a,b] = ksdensity(cmp_pcLactis.corr_coef);
h = area(b,a);
h.EdgeColor = 'k';
h.LineWidth = 1;
h.FaceColor = 'k';
h.FaceAlpha = 0.6;
xlabel('Pearson correlation coefficient','FontSize',14,'FontName','Helvetica');
ylabel('Density','FontSize',14,'FontName','Helvetica');

number = ['N = ',num2str(length(cmp_pcLactis.corr_coef))];
text(-0.9,9,number,'FontSize',15,'FontName','Helvetica','Color','k');

xlim([-1 1]);
set(gcf,'position',[0 0 250 230]);
set(gca,'position',[0.14 0.15 0.85 0.75]);


figure();
hold on;
box on;
[f,x] = ecdf(cmp_pcLactis.corr_coef);
plot(x,f,'LineWidth',1);
clear x f;
set(gca,'FontSize',12,'FontName','Helvetica');

number = ['N = ',num2str(length(cmp_pcLactis.corr_coef))];
text(-0.9,0.9,number,'FontSize',14,'FontName','Helvetica','Color','k');

xlabel('Pearson correlation coefficient','FontSize',14,'FontName','Helvetica');
ylabel('Cumulative frequency','FontSize',14,'FontName','Helvetica');
set(gcf,'position',[500 300 300 300]);
set(gca,'position',[0.2 0.22 0.75 0.75]);

%% Figure distribution of slopes
figure();
hold on;
box on;
[a,b] = ksdensity(sat_factor.k);
h = area(b,a);
h.EdgeColor = 'k';
h.LineWidth = 1;
h.FaceColor = [221,28,119]/256;
h.FaceAlpha = 0.4;
title('Distribution of slope','FontSize',14,'FontName','Helvetica','fontweight','bold');
xlabel('Slope','FontSize',14,'FontName','Helvetica');
ylabel('Density','FontSize',14,'FontName','Helvetica');

number = ['N = ',num2str(length(sat_factor.k))];
text(1.1,2.25,number,'FontSize',15,'FontName','Helvetica','Color','k');

number = ['Median = ',num2str(round(median(sat_factor.k),3))];
text(1.7,2.25,number,'FontSize',15,'FontName','Helvetica','Color',[221,28,119]/256);

set(gcf,'position',[0 500 255 230]);
set(gca,'position',[0.15 0.16 0.8 0.75]);
clear a b h number;


clear;