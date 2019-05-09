%% RNASeqComplement
%   This function returns complementary RNA sequence of a given sequence.
function ComSeq = RNASeqComplement(Seq)
% Input: Seq, a string of a sequence.
% Output: ComSeq, a string of a RNA sequence, complementary to the given
%                 sequence.

ComSeq = '';
n = length(Seq);
for i = 1:n
    Base = Seq(i);
    if strcmpi(Base,'A')
        ComBase = 'U';
    elseif strcmpi(Base,'G')
        ComBase = 'C';
    elseif strcmpi(Base,'C')
        ComBase = 'G';
    elseif strcmpi(Base,'T')
        ComBase = 'A';
    elseif strcmpi(Base,'U')
        ComBase = 'A';
    end
    ComSeq = strcat(ComSeq,ComBase);
end