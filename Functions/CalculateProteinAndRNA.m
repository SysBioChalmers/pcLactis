%% CalculateProteinAndRNA
function [Total_Protein,...
          Total_RNA,...
          Total_rProtein,...
          Total_rRNA,...
          Total_tRNA,...
          Total_mRNA] = CalculateProteinAndRNA(model,sol_full,...
                                               Info_enzyme,...
                                               Info_ribosome,...
                                               Info_tRNA,...
                                               Info_mRNA)

% Calculate f_unmodeled and then weight of the unmodeled protein
MW_unmodeled = 28290.75/1000; %(g/mmol)
rxn_id = strcmp(model.rxns,'R_biomass_dilution');
ump_id = strcmp(model.mets,'unmodeled_protein_biomass_c');
ump_weight = -full(model.S(ump_id,rxn_id)) * MW_unmodeled;

% Growth rate
mu = sol_full(strcmp(model.rxns,'R_biomass_dilution'));

% Ribosome
[~,Total_rProtein,Total_rRNA] = CalculateRibosomeWeight(model,sol_full,Info_ribosome);

% Other proteins (Enzymes)
idx = cellfun(@(x) isequal({'Enzyme dilution'},x),model.subSystems,'UniformOutput',false);
Temp_dil_rxns = model.rxns(cell2mat(idx));
Enzy_dil_rxns = Temp_dil_rxns(~strcmp(Temp_dil_rxns,'R_dilution_ribosome_70S'));
EnzymeID = cellfun(@(x) strcat(x(12:end),'_c'),Enzy_dil_rxns,'UniformOutput',false);
EnzymeWeight = CalculateEnzymeWeight(model,EnzymeID,sol_full,Info_enzyme);

Total_Protein = sum(EnzymeWeight) + Total_rProtein + ump_weight;

% RNA
idx = cellfun(@(x) isequal({'RNA dilution'},x),model.subSystems,'UniformOutput',false);
Temp_dil_rxns = model.rxns(cell2mat(idx));

Total_tRNA = 0;
Total_mRNA = 0;

for i = 1:length(Temp_dil_rxns)
    rxn_id = Temp_dil_rxns{i};
    flux = sol_full(strcmp(model.rxns,rxn_id));

	if contains(rxn_id,'R_dilution_generic_tRNA_')
        trna_name = strrep(rxn_id,'R_dilution_generic_','');
        ge_trna_list = model.rxns(contains(model.rxns,'R_Generic_RNA_'));
        trna_list = ge_trna_list(contains(ge_trna_list,trna_name));
        trna_list = cellfun(@(x) strrep(x,'R_Generic_RNA_',''),...
                                 trna_list,'UniformOutput',false);
        trna_list = cellfun(@(x) x(1:min(strfind(x,trna_name))-2),...
                                 trna_list,'UniformOutput',false);
        MW_list = Info_tRNA.MW(ismember(Info_tRNA.ID,trna_list));
        MW = mean(MW_list)/1000;

        Total_tRNA = Total_tRNA + flux*MW/mu;
        
    else
        tu_name = strrep(rxn_id,'R_dilution_mRNA_','');
        MW = Info_mRNA.MW(strcmp(Info_mRNA.ID,tu_name))/1000;
        
        Total_mRNA = Total_mRNA + flux*MW/mu;       
	end
end

Total_RNA = Total_tRNA + Total_mRNA + Total_rRNA;