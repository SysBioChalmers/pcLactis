%% CalculateRibosomeWeight 
function [Ribosome_weight,...
          rProtein_weight,...
          rRNA_weight]= CalculateRibosomeWeight(model,sol_full,Info_ribosome)
% The unit is g/gCDW.

mu = sol_full(strcmp(model.rxns,'R_biomass_dilution'));

ribo_dil_flux = sol_full(strcmp(model.rxns,'R_dilution_ribosome_70S'));

mol_conc = ribo_dil_flux/mu; %mmol/gCDW

MW_ribosome = Info_ribosome.MW_protein + Info_ribosome.MW_RNA; %g/mol
Ribosome_weight = mol_conc * MW_ribosome / 1000; %g/gCDW

MW_rProtein = Info_ribosome.MW_protein; %g/mol
rProtein_weight = mol_conc * MW_rProtein / 1000; %g/gCDW

MW_rRNA = Info_ribosome.MW_RNA; %g/mol
rRNA_weight = mol_conc * MW_rRNA / 1000; %g/gCDW