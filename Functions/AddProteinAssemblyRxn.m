%% AddProteinAssemblyRxn
function [E_Protein_assembly,...
          Functional_protein] = AddProteinAssemblyRxn(Gene_list)

% Add protein assembly reaction structure.
E_Protein_assembly = MatrixGeneration;

% Import collected protein stoichiometry data.
fid = fopen('protein_stoichiometry_20181030.txt','r');
data = textscan(fid,'%s %f');
fclose(fid);

Protein_stoichiometry(:,1) = data{1};
Protein_stoichiometry(:,2) = num2cell(data{2});

% Main part.
protein_gene_list = Gene_list.total_protein_gene;

Functional_protein = {};

for i = 1:length(protein_gene_list)
    if ismember(protein_gene_list(i,1),Protein_stoichiometry(:,1))%if protein stoichiometry assigned
        index = strcmp(protein_gene_list(i,1),Protein_stoichiometry(:,1));
        stoich = Protein_stoichiometry{index,2};
        if stoich == 1%monomer
            RxnID = strcat('R_',protein_gene_list{i,1},'_Monomer');
            Sub = strcat(protein_gene_list{i,1},'_c');
            Prod = strcat(protein_gene_list{i,1},'_Monomer_c');
            Compo = {Sub;Prod};
            Coeff = [-1;1];
            Rev = 0;
            Subproc = 'Protein assembly';
            Catalyst = '';

            E_Protein_assembly = MatrixOneReactionAdding(E_Protein_assembly,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
            Functional_protein = cat(1,Functional_protein,{Sub(1:end-2),Prod});
        else%2,3,4...mer
            RxnID = strcat('R_',protein_gene_list{i,1},'_',num2str(stoich),'mer');
            Sub = strcat(protein_gene_list{i,1},'_c');
            Prod = strcat(protein_gene_list{i,1},'_',num2str(stoich),'mer_c');            
            Compo = {Sub;Prod};
            Coeff = [-stoich;1];
            Rev = 0;
            Subproc = 'Protein assembly';
            Catalyst = '';

            E_Protein_assembly = MatrixOneReactionAdding(E_Protein_assembly,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
            Functional_protein = cat(1,Functional_protein,{Sub(1:end-2),Prod});
        end
    else%if protein stoichiometry not assigned, then assumed to be 1.
        stoich = 1;
        RxnID = strcat('R_',protein_gene_list{i,1},'_assumed_Monomer');
        Sub = strcat(protein_gene_list{i,1},'_c');
        Prod = strcat(protein_gene_list{i,1},'_assumed_Monomer_c');
        Compo = {Sub;Prod};
        Coeff = [-stoich;1];
        Rev = 0;
        Subproc = 'Protein assembly';
        Catalyst = '';

        E_Protein_assembly = MatrixOneReactionAdding(E_Protein_assembly,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
        Functional_protein = cat(1,Functional_protein,{Sub(1:end-2),Prod});
    end
end
