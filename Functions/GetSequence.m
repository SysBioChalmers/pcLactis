%% GetSequence
%   This function will return sequence of a given gene number.
function sequence = GetSequence(gene,genelist,sequencelist)
% Input: gene, a string of gene ID.
%        genelist, a list of genes.
%        sequencelist, a list of corresponding sequences (AA or NA).
% Output: sequence, the sequence of the given gene, either AA or NA.
index = min(find(strcmp(gene,genelist)));
sequence = cell2mat(sequencelist(index,1));