%% ChangeUnmodeledProteinCluster
%   This function will change the proportion of the unmodeled protein.
function [model,f] = ChangeUnmodeledProteinCluster(model,f_unmodeled)

% Calculate total proportion of proteins (excluding the unmodel protein)
% modeled in the model.
f = 0.46*(1-f_unmodeled); %(g/gCDW)
% f represents the proportion (g/gCDW) of the sum of proteins (without the
% unmodeled protein) considered in the model, which can be used to
% constrain the total fluxes through all dilution reactions of proteins.

% Calculate the stoichiometric coefficient of unmodeled protein in the
% biomass dilution reaction.
MW_unmodeled = 28290.75; %(g/mol)
S_unmodeled = 0.46*f_unmodeled/MW_unmodeled*1000;

rxn_id = strcmp(model.rxns,'R_biomass_dilution');
ump_id = strcmp(model.mets,'unmodeled_protein_biomass_c');
model.S(ump_id,rxn_id) = -S_unmodeled;