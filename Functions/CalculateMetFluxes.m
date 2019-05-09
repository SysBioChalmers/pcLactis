%% CalculateMetFluxes
function met_fluxes = CalculateMetFluxes(model,sol,metid)
% Input: a model, a column of flux distribution and met id.
% Output: a column of atp flux distribution.

met_fluxes = zeros(length(model.rxns),1);

atp_idx = strcmp(model.mets,metid);

for i = 1:length(model.rxns)
    i_flux = sol(i);
    i_coef = model.S(atp_idx,i);
    met_fluxes(i,1) = i_coef * i_flux;
end
 