%% CompareSequence
%   This function checks if two or more genes have identical sequence. All
%   the genes will be compared with the first gene in a given gene list.
function [result, vector] = CompareSequence(genes,Seq_gene_NA_total)
% Input: genes, genes to be compared.
%        Seq_gene_NA_total, a structure containing two columns, i.e., gene and
%                        sequence.
% Output: result, a logical result. Logical 1 means all the genes have
%                 the same sequence. Logical 0 means at least one gene has
%                 different sequence compared to the first gene.
%         vector, a logical vector corresponding to the gene list.
n = length(genes);
vector = zeros(n,1);
seq_1st = GetSequence(genes{1},Seq_gene_NA_total.gene,Seq_gene_NA_total.sequence);
for i = 1:n
    seq = GetSequence(genes{i},Seq_gene_NA_total.gene,Seq_gene_NA_total.sequence);
    vector(i,1) = strcmp(seq_1st,seq);
end
result = all(vector);
