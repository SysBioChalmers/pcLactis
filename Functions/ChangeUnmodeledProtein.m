%% ChangeUnmodeledProtein
%   This function will change the proportion of the unmodeled protein.
function [model,f] = ChangeUnmodeledProtein(model,f_unmodeled)

% Calculate total proportion of RNA and proteins (excluding the unmodel
% protein) modeled in the ME model.
f = 0.46*(1-f_unmodeled)+0.1; %(g/gCDW)

% Calculate the stoichiometric coefficient of unmodeled protein in the
% biomass dilution reaction.
MW_unmodeled = 28290.75; %(g/mol)
S_unmodeled = 0.46*f_unmodeled/MW_unmodeled*1000;

rxn_id = strcmp(model.rxns,'R_biomass_dilution');
ump_id = strcmp(model.mets,'unmodeled_protein_biomass_c');
model.S(ump_id,rxn_id) = -S_unmodeled;