%% AddGenericRNArenamingRxn
function [E_Generic_RNA_renaming,...
          E_RNA_dilution] = AddGenericRNArenamingRxn(FUtRNA,...
                                                     Seq_gene_NA_total,...
                                                     E_RNA_dilution)

% Generate matrix.
E_Generic_RNA_renaming = MatrixGeneration;

% tRNA generic renaming
for i = 1:length(FUtRNA)
    tRNA = FUtRNA{i};
    tRNA_cell = FUtRNA(i);
    pos = max(strfind(tRNA,'_tRNA'));
    genericname = tRNA(pos+1:end-2);

	E_Generic_RNA_renaming = AddRenamingRxn(E_Generic_RNA_renaming,...
                                            tRNA_cell,genericname);
end

% Generate a list of generic tRNA.
index = E_Generic_RNA_renaming.CoeffList == 1;
Generic_RNA = E_Generic_RNA_renaming.CompoList(index);
index_tRNA = ~cellfun('isempty',(cellfun(@(x) strfind(x,'tRNA'),...
                  Generic_RNA,'UniformOutput',false)));
Generic_tRNA = unique(Generic_RNA(index_tRNA));

% tRNA dilution.
for i = 1:length(Generic_tRNA)
    name = Generic_tRNA{i};
    RxnID = strcat('R_dilution_generic_',name(1:end-2));
    Compo = {name};
    Coeff = -1;
    Rev = 0;
    Subproc = 'RNA dilution';
    Catalyst = '';

    E_RNA_dilution = MatrixOneReactionAdding(E_RNA_dilution,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);    
end

% rRNA generic renaming
%All genes that coding the same rRNA have identical sequence, so generic
%rRNA renaming can be done before rRNA modification. However, genes coding
%the same tRNA with reading the same codon may have different sequences, so
%tRNA renaming should be done after tRNA modification.

% Add rRNA.
Gene_list.E_rRNA_16S_gene = {'llmg_rRNA_1';'llmg_rRNA_48a';'llmg_rRNA_4a';
                           'llmg_rRNA_5a';'llmg_rRNA_6b';'llmg_rRNA_7a'};
Gene_list.E_rRNA_23S_gene = {'llmg_rRNA_2';'llmg_rRNA_48c';'llmg_rRNA_4';
                           'llmg_rRNA_5';'llmg_rRNA_6';'llmg_rRNA_7'};
Gene_list.E_rRNA_5S_gene = {'llmg_tRNA_06a';'llmg_rRNA_2a';'llmg_rRNA_48b';
                          'llmg_2426a';'llmg_2466a';'llmg_rRNA_6a';
                          'llmg_rRNA_60a'};

% Check if genes coding the same rRNA have identical sequence.
[res16,~] = CompareSequence(Gene_list.E_rRNA_16S_gene,Seq_gene_NA_total);
[res23,~] = CompareSequence(Gene_list.E_rRNA_23S_gene,Seq_gene_NA_total);
if res16 == 0
    error('Sequences in 16S rRNA are different.');
end
if res23 == 0
    error('Sequences in 23S rRNA are different.');
end

% Considering the fact that all 16S rRNA genes have identical sequence and
% all 23S rRNA also have identical sequence, generic reactions are made for
% 16S rRNA and 23S rRNA before rRNA modification. This decreases the number
% of reactions. 5S rRNA genes have different sequences but no modification 
% reactions will be formulated for them, so generic reactions can also made
% for them in this stage.
%   Name component ID in the model.
rRNA_16S = cellfun(@(x) strcat(x,'_unmodified_c'),...
                   Gene_list.E_rRNA_16S_gene,'UniformOutput',false);
rRNA_23S = cellfun(@(x) strcat(x,'_unmodified_c'),...
                   Gene_list.E_rRNA_23S_gene,'UniformOutput',false);
rRNA_5S = cellfun(@(x) strcat(x,'_unmodified_c'),...
                   Gene_list.E_rRNA_5S_gene,'UniformOutput',false);
%   Generate matrix for generic rRNA renaming.
E_Generic_RNA_renaming = AddRenamingRxn(E_Generic_RNA_renaming,...
                           rRNA_16S,'rRNA_16S_unmodified');
E_Generic_RNA_renaming = AddRenamingRxn(E_Generic_RNA_renaming,...
                           rRNA_23S,'rRNA_23S_unmodified');
E_Generic_RNA_renaming = AddRenamingRxn(E_Generic_RNA_renaming,...
                           rRNA_5S,'rRNA_5S_unmodified');
clear rRNA_16S rRNA_23S rRNA_5S;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_Generic_RNA_renaming = AddRenamingRxn(E_Generic_RNA_renaming,...
                                                 genelist,...
                                                 genericname)

for i = 1:length(genelist)
    gene = genelist{i};
    gene_process = gene(1:end-2);
    RxnID = strcat('R_Generic_RNA_',gene_process,'_to_',genericname);
    Compo = {genelist{i};strcat(genericname,'_c')};
    Coeff = [-1;1];
    Rev = 0;
    Subproc = 'Generic RNA renaming';
    Catalyst = '';

    E_Generic_RNA_renaming = MatrixOneReactionAdding(E_Generic_RNA_renaming,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
end
end
