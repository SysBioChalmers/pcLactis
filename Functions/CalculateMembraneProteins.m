%% CalculateMembraneProteins
function f_memprot = CalculateMembraneProteins(model,sol_full,Info_enzyme)
% Input: a column of flux distribution.
% Output: proportion of membrane proteins (g/g biomass).

mu = sol_full(contains(model.rxns,'R_biomass_dilution'));

Transporters = {};
idx = cellfun(@(x) isequal({'Metabolism'},x),model.subSystems,'UniformOutput',false);
M_rxns_idx = find(cell2mat(idx));
for i = 1:length(M_rxns_idx)
    rxnidx = M_rxns_idx(i);
	metlist = model.mets(model.S(:,rxnidx) ~= 0);
	if find(cell2mat(cellfun(@(x) strcmp(x(length(x)-1:end),'_e'),metlist,'UniformOutput',false)),1) > 0
        enzyme_id = model.grRules{rxnidx};
        if ~isempty(enzyme_id)
            Transporters = [Transporters;{enzyme_id}];
        end
	end
end

f_memprot = 0;
for i = 1:length(Transporters)
	comp_name = Transporters{i};
	rxn_id = strcat('R_dilution_',comp_name(1:length(comp_name)-2));
	flux = sol_full(strcmp(model.rxns,rxn_id));

	MW = Info_enzyme.MW(contains(Info_enzyme.ID,comp_name));
	coeff = MW/1000;
    f_i = coeff*flux/mu;
    f_memprot = f_memprot + f_i;
end