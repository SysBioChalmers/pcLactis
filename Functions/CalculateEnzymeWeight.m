%% CalculateEnzymeWeight 
function EnzymeWeight = CalculateEnzymeWeight(model,EnzymeID,...
                                              sol_full,Info_enzyme)
% The unit is g/gCDW. EnzymeID should be in cell structure.

mu = sol_full(strcmp(model.rxns,'R_biomass_dilution'));

EnzymeWeight = [];

for i = 1:length(EnzymeID)
    id = EnzymeID{i};
    dil_rxn_id = strcat('R_dilution_',id(1:end-2));
    dil_flux = sol_full(strcmp(model.rxns,dil_rxn_id));
    mol_conc = dil_flux/mu; %mmol/gCDW
    MW = Info_enzyme.MW(strcmp(Info_enzyme.ID,id)); %g/mol
    weight = mol_conc * MW / 1000; %g/gCDW
    EnzymeWeight = [EnzymeWeight;weight];
end
