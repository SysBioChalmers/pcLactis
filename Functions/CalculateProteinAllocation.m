%% CalculateProteinAllocation
function [pathway, absvalue] = CalculateProteinAllocation(model,sol_full,Info_enzyme,excel_input)

pathway = cell(0,1);
absvalue = zeros(0,1);

[~, n] = size(excel_input);
for i = 1:n
    pathway(i,1) = excel_input(1,i);
    enzymes_tmp = excel_input(2:end,i);
    enzymes = enzymes_tmp(contains(enzymes_tmp,'_Enzyme_c'));
    absvalue(i,1) = sum(CalculateEnzymeWeight(model,enzymes,sol_full,Info_enzyme));
end
