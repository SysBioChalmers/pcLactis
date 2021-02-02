% Remove blocked enzymes

load('M_Model.mat');
model = M_Model;

blockedRxns = findBlockedReaction(model);
blockedRxns = blockedRxns';
blockedRxns = cellfun(@(x) x(3:end),blockedRxns,'UniformOutput',false);
blockedRxns = cellfun(@(x) strcat(x,'_'),blockedRxns,'UniformOutput',false);

load('kcat_M.mat');
m_enzyme = kcat_M.gpr;
[~, ~, k_raw_e] = xlsread('kcat_E.xlsx','E');
e_enzyme = k_raw_e(2:end,1);

param_list = [m_enzyme;e_enzyme];

idx = contains(param_list,blockedRxns);

param_list = param_list(~idx);

cd Cluster/;
save('param_list.mat','param_list');
cd ../;
