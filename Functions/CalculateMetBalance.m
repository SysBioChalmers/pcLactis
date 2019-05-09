%% CalculateMetBalance
function netflux = CalculateMetBalance(model,rxnlist,sol,metid)
% Input: a model, a list of rxn ids with fluxes and met id.
% Output: net flux of the metabolite in the rxn list.

coeff = model.S(strcmp(model.mets,metid),ismember(model.rxns,rxnlist));
netflux = coeff * sol;
