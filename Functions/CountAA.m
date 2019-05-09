%% CountAA
%   This function will count the number of amino acid in a given sequence.
function sum = CountAA(AAsequence)
% Input: AAsequence, a string of sequence made up of AA abbr.
% Output: sum, a structure containing three fields, i.e., AA abbr, met ID
%         and number

cd Data;
AA_ID = importdata('AA_ID.txt');
cd ../;
AA_ID = split(AA_ID);
sum = struct();
sum.metid = AA_ID(:,2);
sum.abbr = AA_ID(:,3);
sum.num = [];

for i = 1:length(sum.abbr)
    AA_abbr = sum.abbr{i};
    sum.num(i,1) = length(strfind(AAsequence,AA_abbr));
end

