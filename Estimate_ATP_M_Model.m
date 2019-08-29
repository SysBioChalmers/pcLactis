%% Basic setting
load('M_Model.mat');
model = M_Model;

% Block uptake direction of all the exchange reactions.
idx = cellfun(@(x) contains(x,'R_M_EX_'),model.rxns,'UniformOutput',false);
M_exchange_reactions = model.rxns(cell2mat(idx));
model = changeRxnBounds(model,M_exchange_reactions,0,'l');
clear idx M_exchange_reactions;

% Block some reactions in M part.
model = changeRxnBounds(model,'R_M_biomass_LLA',0,'b');
% model = changeRxnBounds(model,'R_M_MGt2pp',0,'l');%block infinite h[e]

% Set NGAM.
model=changeRxnBounds(model,'R_M_ATPM',0,'l');
model=changeRxnBounds(model,'R_M_ATPM',1000,'u');

% Set objective funtion.
model = changeObjective(model,'R_M_ATPM',1);

%% Loop for calculating max ATP of each experiment
% Set bounds for some exchange reactions.
[~, ~, exchange_raw] = xlsread('Exchange_reaction_setting.xlsx','Exp_bounds_GAM_NGAM');
Exchange_reactions = exchange_raw(2:end,1);
col_header = exchange_raw(1,1:21);
max_atp = zeros((length(col_header)-1)/2,2);
for i = 2:2:length(col_header)
    header = col_header{i};
    mu = str2double(header(4:7));
    LB = cell2mat(exchange_raw(2:end,i));
    UB = cell2mat(exchange_raw(2:end,i+1));
    model = changeRxnBounds(model,Exchange_reactions,LB,'l');
    model = changeRxnBounds(model,Exchange_reactions,UB,'u');
    model = changeRxnBounds(model,'R_M_biomass_LLA_atpm',mu,'b');
    sol = optimizeCbModel(model,'max','one');
    max_atp(i/2,:) = [mu,sol.f];
end