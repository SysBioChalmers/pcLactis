%% AddUnmodeledProtein
%   This function will add a unmodel protein into the model. Here is used
%   for a special protein, i.e., an unmodeled protein which represents the
%   rest of the proteins not modeled in the model. And the synthesis of
%   this protein also accounts for the RNA usuage.
function [Gene_list,...
          TU,...
          Seq_gene_NA_total,...
          Seq_CDS_AA_total,...
          E_Protein_assembly,...
          Functional_protein,...
          E_Enzyme_formation,...
          E_Protein_degradation,...
          E_Enzyme_dilution] = AddUnmodeledProtein(Gene_list,...
                                                   TU,...
                                                   Seq_gene_NA_total,...
                                                   Seq_CDS_AA_total,...
                                                   E_Protein_assembly,...
                                                   Functional_protein,...
                                                   E_Enzyme_formation,...
                                                   E_Protein_degradation,...
                                                   E_Enzyme_dilution)

% Import seq info of the unmodel protein.
[~, ~, seq_raw] = xlsread('Unmodeled_protein.xlsx','Seq_info_unmodeled');

AA_sequence = 'M';
NA_sequence = 'ATG';

for i = 2:length(seq_raw)-1
    AA_Abbr = seq_raw{i,1};
    num_AA = seq_raw{i,2};
    Codon = seq_raw{i,3};
    AA_sequence = strcat(AA_sequence,repmat(AA_Abbr,1,num_AA));
    NA_sequence = strcat(NA_sequence,repmat(Codon,1,num_AA));
end
NA_sequence = strcat(NA_sequence,seq_raw{end,3});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input basic info.
gene_name = 'unmodeled';
complex_name = 'unmodeled_protein';
TU_name = 'TUunmodeled_1';
NA_seq = NA_sequence;
AA_seq = AA_sequence;
TU_seq = NA_seq;
protein_stoichiometry = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add the gene name to the list Gene_list.total_protein_gene.
Gene_list.total_protein_gene = cat(1,Gene_list.total_protein_gene,{gene_name});

% Add the TU info to the structure TU.
TU.TUID = cat(1,TU.TUID,{TU_name});
TU.gene = cat(1,TU.gene,{gene_name});
TU.sequence = cat(1,TU.sequence,{TU_seq});

% Add the TU name to the list Gene_list.TU.
Gene_list.TU = cat(1,Gene_list.TU,{TU_name});

% Add the NA seq to the structure Seq_gene_NA_total.
Seq_gene_NA_total.gene = cat(1,Seq_gene_NA_total.gene,{gene_name});
Seq_gene_NA_total.sequence = cat(1,Seq_gene_NA_total.sequence,{NA_seq});

% Add the AA seq to the structure Seq_gene_NA_total.
Seq_CDS_AA_total.gene = cat(1,Seq_CDS_AA_total.gene,{gene_name});
Seq_CDS_AA_total.sequence = cat(1,Seq_CDS_AA_total.sequence,{AA_seq});

% Add protein assembly reaction.
if protein_stoichiometry == 1%monomer
	RxnID = strcat('R_',gene_name,'_Monomer');
	Sub = strcat(gene_name,'_c');
	Prod = strcat(gene_name,'_Monomer_c');
	Compo = {Sub;Prod};
	Coeff = [-1;1];
	Rev = 0;
	Subproc = 'Protein assembly';
	Catalyst = '';
    E_Protein_assembly = MatrixOneReactionAdding(E_Protein_assembly,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
	Functional_protein = cat(1,Functional_protein,{Sub(1:end-2),Prod});
else%2,3,4...mer
	RxnID = strcat('R_',gene_name,'_',num2str(protein_stoichiometry),'mer');
	Sub = strcat(gene_name,'_c');
	Prod = strcat(gene_name,'_',num2str(protein_stoichiometry),'mer_c');            
	Compo = {Sub;Prod};
	Coeff = [-protein_stoichiometry;1];
	Rev = 0;
	Subproc = 'Protein assembly';
	Catalyst = '';
    E_Protein_assembly = MatrixOneReactionAdding(E_Protein_assembly,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
	Functional_protein = cat(1,Functional_protein,{Sub(1:end-2),Prod});
end

% Add complex formation reaction.
RxnID = strcat('R_',complex_name);
Compo = {Prod;strcat(complex_name,'_c')};
Coeff = [-1;1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                 RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

% Add protein degradation reaction.    
RxnID = strcat('R_',complex_name,'_degradation');
Compo = {strcat(complex_name,'_c');strcat(Prod(1:end-2),'_degradation_c')};
Coeff = [-1;1];
Rev = 0;
Subproc = 'Protein degradation';
Catalyst = '';
E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

% Add dilution reaction.
% A dilution reaction should also be added for the protein, but for the
% unmodeled protein, the dilution reaction will be integrated into the
% biomass equation. 
RxnID = strcat('R_dilute_',complex_name,'_to_biomass');
Compo = {strcat(complex_name,'_c');strcat(complex_name,'_biomass_c')};
Coeff = [-1;1];
Rev = 0;
Subproc = '';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);