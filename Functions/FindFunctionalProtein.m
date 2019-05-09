%% FindFunctionalProtein
%   This function returns the functional protein ID of the gene ID.
function Functional_protein_ID = FindFunctionalProtein(gene_ID,Functional_protein_list)
% Input: gene_ID, a gene will be searched.
%        Functional_protein_list, a list containing gene&functional protein
%                                 relationship.
% Output: Functional_protein_ID, the corresponding functional protein ID.
if ismember(gene_ID,Functional_protein_list(:,1))
    index = strcmp(gene_ID,Functional_protein_list(:,1));
    Functional_protein_ID = cell2mat(Functional_protein_list(index,2));
else
    error('The gene you are retrieving is not in the functional protein list.');
end
end