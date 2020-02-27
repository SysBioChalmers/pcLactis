%% CalculateAAtRNAFlux
function flux = CalculateAAtRNAFlux(model,sol,aaid)
% Calculate tRNA charging flux for a given amino acid.

idx_tRNAcharging = contains(model.rxns,'_charging_tRNA_');
idx_withGPR = contains(model.grRules,'tRNA_Synthetase_');
idx = idx_tRNAcharging & idx_withGPR;

rxnlist = model.rxns(idx);

AArxn = rxnlist(contains(lower(rxnlist),strcat(aaid,'_')));
flux = sum(sol(ismember(model.rxns,AArxn),:),1);

