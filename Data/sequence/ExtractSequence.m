%% ExtractSequence
%   This function processes sequence information from NCBI, and imports as
%   a .mat structure.
function seq = ExtractSequence(file)
% Input: file, a text file downloaded from NCBI.
% Output: a structure containing geneID and corresponding sequence.
raw = importdata(file);
seq = struct();
seq.gene = cell(0,1);
seq.seq = cell(0,1);
n = find(startsWith(raw,'>'));
for i = 1:length(n)
    x = length(find(cellfun('isempty',seq.gene) == 0));
    pos = n(i);
    gene = raw(pos);
    gene = extractBetween(gene,'[locus_tag=',']');
    seq.gene(x+1,1) = gene;
end
m = [n;length(raw)+1];
for i = 1:length(m)-1
    x = length(find(cellfun('isempty',seq.seq) == 0));
    sequence_rows = transpose([m(i)+1:m(i+1)-1]);
    sequence_in_row = raw(sequence_rows,1);
    sequence = sequence_in_row{1};
    for j = 1:length(sequence_in_row)-1
        sequence = strcat(sequence,sequence_in_row{j+1});
    end
    seq.seq(x+1,1) = {sequence};
end
