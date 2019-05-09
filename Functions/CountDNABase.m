%% CountDNABase
%   This function will count the number of DNA base in a given sequence.
function bases = CountDNABase(DNAsequence)
% Input: DNAsequence, a string of sequence made up of A, G, C and T.
% Output: sum, a vector containing the number of each base in the sequence.
%         sum = [A G C T]
%         A, the number of adenine in the sequence.
%         G, the number of guanine in the sequence.
%         C, the number of cytosine in the sequence.
%         T, the number of thymine in the sequence.
A = length(strfind(DNAsequence,'A'));
G = length(strfind(DNAsequence,'G'));
C = length(strfind(DNAsequence,'C'));
T = length(strfind(DNAsequence,'T'));
bases = [A G C T];
