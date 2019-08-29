%% ProcessSequence
%   This function processes sequences downloaded from databases.
function [Seq_CDS_AA_total, Seq_gene_NA_total, Gene_list] = ProcessSequence(Gene_list)

Chr_CDS_AA = ExtractSequence('Chromosome_CDS_protein_sequence.txt');
Chr_gene_NA = ExtractSequence('Chromosome_gene_nucleotide_sequence.txt');
Seq_CDS_AA = struct();
Seq_CDS_AA.gene = Chr_CDS_AA.gene;
Seq_CDS_AA.sequence = Chr_CDS_AA.seq; 
Seq_gene_NA = struct();
Seq_gene_NA.gene = Chr_gene_NA.gene;
Seq_gene_NA.sequence = Chr_gene_NA.seq;

% Count the number of total protein genes based on TUs.
%That is to remove non-CDS genes from Gene_list.total_gene_in_all_TUs.

non_CDS_genes = Gene_list.E_RNA_gene;      

Gene_list.total_protein_gene_in_all_TUs = {};
for i = 1:length(Gene_list.total_gene_in_all_TUs)
    gene_name = Gene_list.total_gene_in_all_TUs{i};
    if ~ismember(gene_name,non_CDS_genes)
        Gene_list.total_protein_gene_in_all_TUs = cat(1,...
                                  Gene_list.total_protein_gene_in_all_TUs,...
                                  {gene_name});
    end
end

% Extract sequence for genes in the ME model.
% ID conversion will be done below as gene ID used in NCBI is different
% (which is newer) from the ID used in the model. The ID relationship was
% downloaded from BioCyc. For some missing ID relationships we assigned
% manually.
[~, ~, IDconvertion] = xlsread('Locus_ID_convertion.xlsx');
newID = IDconvertion(:,1);
oldID = IDconvertion(:,2);

Seq_CDS_AA_total = struct();
Seq_CDS_AA_total.gene = Gene_list.total_protein_gene_in_all_TUs;
Seq_CDS_AA_total.sequence = cell(0,1);
for i = 1:length(Seq_CDS_AA_total.gene)
    x = length(find(cellfun('isempty',Seq_CDS_AA_total.sequence) == 0));
    ID_ME_old = Seq_CDS_AA_total.gene{i};
	if find(strcmp(oldID,ID_ME_old)) > 0
        index_1 = strcmp(oldID,ID_ME_old);
        ID_ME_new = newID{index_1};
        if ismember(ID_ME_new,Seq_CDS_AA.gene)
            index_2 = strcmp(Seq_CDS_AA.gene,ID_ME_new);
            Seq_CDS_AA_total.sequence(x+1,1) = Seq_CDS_AA.sequence(index_2,1);
        else
            error(['The gene ',ID_ME_new,' is not in Seq_CDS_AA.gene.']);
        end
    else
        error(['The gene ',ID_ME_old,' is not in Locus_ID_convertion.xlsx.']);
	end
end
Seq_gene_NA_total = struct();
Seq_gene_NA_total.gene = Gene_list.total_gene_in_all_TUs;
Seq_gene_NA_total.sequence = cell(0,1);
for i = 1:length(Seq_gene_NA_total.gene)
    x = length(find(cellfun('isempty',Seq_gene_NA_total.sequence) == 0));
    ID_ME_old = Seq_gene_NA_total.gene{i};
	if find(strcmp(oldID,ID_ME_old)) > 0
        index_1 = strcmp(oldID,ID_ME_old);
        ID_ME_new = newID{index_1};
        if ismember(ID_ME_new,Seq_gene_NA.gene)
            index_2 = strcmp(Seq_gene_NA.gene,ID_ME_new);
            Seq_gene_NA_total.sequence(x+1,1) = Seq_gene_NA.sequence(index_2,1);
        else
            error(['The gene ',ID_ME_new,' is not in Seq_gene_NA.gene.']);
        end
	else
        error(['The gene ',ID_ME_old,' is not in Locus_ID_convertion.xlsx.']);
	end
end