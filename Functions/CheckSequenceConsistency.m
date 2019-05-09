%% CheckSequenceConsistency
%   This function checks if all the genes in the model have intact
%   sequences involved in transcription units (TUs).
function [result, vector] = CheckSequenceConsistency(genelist,TU)
% Input: genelist, a structure containing two columns, one is gene and the
%                  other is sequence. The sequence should be nucleic acid 
%                  sequence.
%        TU, a structure containing three columns, i.e., TUID, gene and
%            sequence.
% Output: result, a logical result. Logical 1 means all the genes have
%                 intact sequences in TUs. Logical 0 means at least one 
%                 gene does not have an intact sequence in any TUs.
%         vector, a logical vector corresponding to the gene list.
n = length(genelist.gene);
vector = zeros(n,1);
for i = 1:n
    gene = genelist.gene{i};
    gene_seq = genelist.sequence{i};
    TU_index = find(strcmp(TU.gene,gene));
    m = length(TU_index);
    if m == 1
        TU_seq = TU.sequence{TU_index};
        vector(i,1) = ~isempty(strfind(TU_seq,gene_seq));
    elseif m > 1
        for j = 1:m
            subvector = zeros(m,1);
            subindex = TU_index(j);
            TU_seq = TU.sequence{subindex};
            subvector(j,1) = ~isempty(strfind(TU_seq,gene_seq));
        end
        vector(i,1) = all(subvector);
    end
end
result = all(vector);
