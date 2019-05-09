%% ProcessSequence
%   This function processes sequences downloaded from databases.
function [Seq_CDS_AA_total, Seq_gene_NA_total, Gene_list] = ProcessSequence(Gene_list)

cd Data/sequence;
Chr_CDS_AA = ExtractSequence('Chromosome_CDS_protein_sequence.txt');
Chr_gene_NA = ExtractSequence('Chromosome_gene_nucleotide_sequence.txt');
pLP712_CDS_AA = ExtractSequence('pLP712_CDS_protein_sequence.txt');
pLP712_gene_NA = ExtractSequence('pLP712_gene_nucleotide_sequence.txt');
pNZ712_CDS_AA = ExtractSequence('pNZ712_CDS_protein_sequence.txt');
pNZ712_gene_NA = ExtractSequence('pNZ712_gene_nucleotide_sequence.txt');
Seq_CDS_AA = struct();
Seq_CDS_AA.gene = [Chr_CDS_AA.gene;pLP712_CDS_AA.gene;pNZ712_CDS_AA.gene];
Seq_CDS_AA.sequence = [Chr_CDS_AA.seq;pLP712_CDS_AA.seq;pNZ712_CDS_AA.seq]; 
Seq_gene_NA = struct();
Seq_gene_NA.gene = [Chr_gene_NA.gene;pLP712_gene_NA.gene;pNZ712_gene_NA.gene];
Seq_gene_NA.sequence = [Chr_gene_NA.seq;pLP712_gene_NA.seq;pNZ712_gene_NA.seq];

%It was found that the gene 'pNZ712_25' was assigned with a wrong NA
%sequence, i.e., the number of NAs is not the multiple of 3. The following
%code is used to fix this.
[~, index] = ismember('pNZ712_25',Seq_gene_NA.gene);
seq_pNZ712_25 = Seq_gene_NA.sequence{index};
n = length(seq_pNZ712_25)/3;
if n ~= fix(n)
    Seq_gene_NA.sequence{index} = seq_pNZ712_25(1:end-1);
end

% Count the number of total protein genes based on TUs.
%That is to remove non-CDS genes from Gene_list.total_gene_in_all_TUs. The
%non-CDS genes consist of E_RNA_gene and non-CDS genes in plasmid.
non_CDS_genes_plasmid = {};
pl_gene = unique(cat(1,pLP712_gene_NA.gene,pNZ712_gene_NA.gene));
pl_CDS_gene = unique(cat(1,pLP712_CDS_AA.gene,pNZ712_CDS_AA.gene));
for i = 1:length(pl_gene)
    gene_name = pl_gene{i};
    if ~ismember(gene_name,pl_CDS_gene)
        non_CDS_genes_plasmid = cat(1,non_CDS_genes_plasmid,{gene_name});
    end
end

non_CDS_genes = unique(cat(1,non_CDS_genes_plasmid,Gene_list.E_RNA_gene));      

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
    if startsWith(ID_ME_old,'llmg')%no need to convert plasmid gene IDs.
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
    else%for plasmid genes
        pID = Seq_CDS_AA_total.gene{i};
        if ismember(pID,Seq_CDS_AA.gene)
            index_p = strcmp(Seq_CDS_AA.gene,pID);
            Seq_CDS_AA_total.sequence(x+1,1) = Seq_CDS_AA.sequence(index_p,1);
        else
            error(['The gene ',pID,' is not in Seq_CDS_AA.gene.']);
        end
    end
end
Seq_gene_NA_total = struct();
Seq_gene_NA_total.gene = Gene_list.total_gene_in_all_TUs;
Seq_gene_NA_total.sequence = cell(0,1);
for i = 1:length(Seq_gene_NA_total.gene)
    x = length(find(cellfun('isempty',Seq_gene_NA_total.sequence) == 0));
    ID_ME_old = Seq_gene_NA_total.gene{i};
    if startsWith(ID_ME_old,'llmg')%no need to convert plasmid gene IDs.
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
    else%for plasmid genes
        pID = Seq_gene_NA_total.gene{i};
        if ismember(pID,Seq_gene_NA.gene)
            index_p = strcmp(Seq_gene_NA.gene,pID);
            Seq_gene_NA_total.sequence(x+1,1) = Seq_gene_NA.sequence(index_p,1);
        else
            error(['The gene ',pID,' is not in Seq_CDS_AA.gene.']);
        end
    end
end
cd ../../;