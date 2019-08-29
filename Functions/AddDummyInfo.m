%% AddDummyInfo
%   This function will add information for the dummy gene.
function [Gene_list,...
          TU,...
          Seq_gene_NA_total,...
          Seq_CDS_AA_total] = AddDummyInfo(Gene_list,...
                                           TU,...
                                           Seq_gene_NA_total,...
                                           Seq_CDS_AA_total)

% Import seq info of the unmodeled protein.
[~, ~, seq_raw] = xlsread('Dummy_GPR.xlsx','Seq_info_dummy');

AA_sequence = 'M';
NA_sequence = 'ATG';

for i = 2:length(seq_raw)-1
    AA_Abbr = seq_raw{i,1};
    num_AA = seq_raw{i,2};
    Codon = seq_raw{i,3};
    AA_sequence = strcat(AA_sequence,repmat(AA_Abbr,1,num_AA));
    NA_sequence = strcat(NA_sequence,repmat(Codon,1,num_AA));
end
NA_sequence = strcat(NA_sequence,seq_raw{end,3});

% Input basic info.
gene_name = 'dummy';
TU_name = 'TUdummy_1';
NA_seq = NA_sequence;
AA_seq = AA_sequence;
TU_seq = NA_seq;

% Add the gene name to the list Gene_list.total_protein_gene.
if isempty(find(strcmp(Gene_list.total_protein_gene,gene_name), 1))
    Gene_list.total_protein_gene = cat(1,Gene_list.total_protein_gene,{gene_name});
end
% Add the TU info to the structure TU.
if isempty(find(strcmp(TU.TUID,TU_name), 1))
    TU.TUID = cat(1,TU.TUID,{TU_name});
    TU.gene = cat(1,TU.gene,{gene_name});
    TU.sequence = cat(1,TU.sequence,{TU_seq});
end

% Add the TU name to the list Gene_list.TU.
if isempty(find(strcmp(Gene_list.TU,TU_name), 1))
    Gene_list.TU = cat(1,Gene_list.TU,{TU_name});
end

% Add the NA seq to the structure Seq_gene_NA_total.
Seq_gene_NA_total.gene = cat(1,Seq_gene_NA_total.gene,{gene_name});
Seq_gene_NA_total.sequence = cat(1,Seq_gene_NA_total.sequence,{NA_seq});

% Add the AA seq to the structure Seq_gene_NA_total.
Seq_CDS_AA_total.gene = cat(1,Seq_CDS_AA_total.gene,{gene_name});
Seq_CDS_AA_total.sequence = cat(1,Seq_CDS_AA_total.sequence,{AA_seq});


