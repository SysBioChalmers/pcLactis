%% CalculateTranscriptomics 
function AbsTranscripts = CalculateTranscriptomics(model,sol_full)
% The unit of concentration is mmol/gCDW.

[~, txtTUinfo, ~] = xlsread('Transcription_units.xlsx');

TUID = txtTUinfo(2:end,1);
TUgene = txtTUinfo(2:end,2);

idx = cellfun(@(x) isequal({'RNA dilution'},x),model.subSystems,...
              'UniformOutput',false);
rna_rxns = model.rxns(cell2mat(idx));
rna_dil_rxns = rna_rxns(~contains(rna_rxns,'_generic_tRNA_'));

TU_involved = cellfun(@(x) strrep(x,'R_dilution_mRNA_',''),rna_dil_rxns,...
                      'UniformOutput',false);
TUlist = TUID(ismember(TUID,TU_involved));
genelist = TUgene(ismember(TUID,TU_involved));
conclist = zeros(length(genelist),1);

mu = sol_full(strcmp(model.rxns,'R_biomass_dilution'));

for i = 1:length(genelist)
    TU_id = TUlist{i};
    dil_rxn_id = strcat('R_dilution_mRNA_',TU_id);
    flux = sol_full(strcmp(model.rxns,dil_rxn_id));
    conclist(i,1) = flux/mu;
end

AbsTranscripts = struct();
AbsTranscripts.gene_id = unique(genelist);
AbsTranscripts.concentration = zeros(length(AbsTranscripts.gene_id),1);

for i = 1:length(AbsTranscripts.gene_id)
    geneid = AbsTranscripts.gene_id(i);
    AbsTranscripts.concentration(i,1) = sum(conclist(strcmp(geneid,genelist)));
end

