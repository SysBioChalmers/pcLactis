%% AddProteinMaturationRxn
%   This function will generate protein maturation reactions.
function [E_Protein_maturation,...
          Translated_gene,...
          E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = AddProteinMaturationRxn(Translated_gene,...
                                                           Functional_protein,...
                                                           E_Enzyme_formation,...
                                                           E_Enzyme_dilution,...
                                                           E_Protein_degradation)

% Generate matrix.
E_Protein_maturation = MatrixGeneration;

% Check if every protein starts with Met.
for i = 1:length(Translated_gene.AAseq)
    gene = Translated_gene.ID{i};
    seq = Translated_gene.AAseq{i};
    if ~startsWith(seq,'M')
        error(['The protein sequence of ',gene,' does not start with Met.']);
    end    
end

% Add peptide deformylase.
%Def -> Def
RxnID = 'R_Peptide_Deformylase_Enzyme';
Def = FindFunctionalProtein('llmg_0532',Functional_protein);
Product = 'Peptide_Deformylase_Enzyme_c';
Compo = {Def;Product};
Coeff = [-1;1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

%add protein degradation reactions
RxnID = 'R_Peptide_Deformylase_Enzyme_degradation';
Sub = strcat(Def(1:end-2),'_degradation_c');
Compo = {Product;Sub};
Coeff = [-1;1];
Rev = 0;
Subproc = 'Protein degradation';
Catalyst = '';
E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                         
%add dilution reactions
RxnID = 'R_dilution_Peptide_Deformylase_Enzyme';
Compo = {Product};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 
% Add methionine aminopeptidase.
%Map -> Map
RxnID = 'R_Methionine_Aminopeptidase_Enzyme';
Map = FindFunctionalProtein('llmg_0577',Functional_protein);
Product = 'Methionine_Aminopeptidase_Enzyme_c';
Compo = {Map;Product};
Coeff = [-1;1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);                                 

%add protein degradation reactions
RxnID = 'R_Methionine_Aminopeptidase_Enzyme_degradation';
Sub = strcat(Map(1:end-2),'_degradation_c');
Compo = {Product;Sub};
Coeff = [-1;1];
Rev = 0;
Subproc = 'Protein degradation';
Catalyst = '';
E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                         
%add dilution reactions
RxnID = 'R_dilution_Methionine_Aminopeptidase_Enzyme';
Compo = {Product};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 
% N-terminus prediction
% The prediction of N-terminal protein was perform using the TermiNator
% (https://bioweb.i2bc.paris-saclay.fr/terminator3/), stored in the file
% 'N_terminal_prediction.xlsx'. The proteins with likelyhood < 80% will be
% removed in the following codes, and then it will be assumed that their
% mature proteins start with Met, i.e., the pos value will be assumed to
% be 1. This means that no excision will be done.
cd Data;
[numNT, txtNT, ~] = xlsread('N_terminus_prediction.xlsx');
cd sequence;
[~, txtID, ~] = xlsread('Locus_ID_convertion.xlsx');
cd ../../;

ID_new = txtID(:,1);
ID_old = txtID(:,2);
NT_IDnew = txtNT(2:end,1);
NT_lik = numNT;
NT_pos = txtNT(2:end,2);
NT_pos_processed = [];
for i = 1:length(NT_pos)
    pos_raw = NT_pos{i};
    pos_process = str2num(cell2mat((extractBetween(pos_raw,'(',')'))));
    NT_pos_processed = [NT_pos_processed;pos_process];
end

% Collect AAseq of mature proteins.
Translated_gene.AAseq_mature = {};

% Add main reactions.
for i = 1:length(Translated_gene.ID)
    gene_ID = Translated_gene.ID{i};
    pos = NTerminusPrediction(gene_ID,ID_new,ID_old,NT_IDnew,NT_lik,...
                              NT_pos_processed);
    if pos == 1%1st amino acid of the mature protein is Met
        Translated_gene.AAseq_mature{i,1} = Translated_gene.AAseq{i};
        % Remove formyl-group
        sub_peptide = strcat(gene_ID,'_nascent_c');
        prod_peptide = strcat(gene_ID,'_c');
        E_Protein_maturation = AddDefRxn(E_Protein_maturation,...
                                         gene_ID,...
                                         sub_peptide,prod_peptide);
    else%Met should be removed
        if pos == 2
            Translated_gene.AAseq_mature{i,1} = Translated_gene.AAseq{i}(2:end);
            % Remove formyl-group
            sub_peptide = strcat(gene_ID,'_nascent_c');
            prod_peptide = strcat(gene_ID,'_for_M_excision_c');  
            E_Protein_maturation = AddDefRxn(E_Protein_maturation,...
                                             gene_ID,...
                                             sub_peptide,prod_peptide);

            % Remove methionyl-group
            sub_peptide = strcat(gene_ID,'_for_M_excision_c');
            prod_peptide = strcat(gene_ID,'_c');
            E_Protein_maturation = AddMapRxn(E_Protein_maturation,...
                                             gene_ID,...
                                             sub_peptide,prod_peptide);
        else%The first AA of the mature protein is not the second AA of the nascent peptide
            error(['More than the 1st Met should be removed for ',gene_ID,'.']);
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_Protein_maturation = AddDefRxn(E_Protein_maturation,...
                                      gene_ID,...
                                      sub_peptide,prod_peptide)
%peptide-nascent + h2o -> peptide + for
RxnID = strcat('R_deformylase_',gene_ID);
Compo = {sub_peptide;'M_h2o_c';prod_peptide;'M_for_c'};
Coeff = [-1;-1;1;1];
Rev = 0;
Subproc = 'Protein maturation';
Catalyst = 'Peptide_Deformylase_Enzyme_c';

E_Protein_maturation = MatrixOneReactionAdding(E_Protein_maturation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_Protein_maturation = AddMapRxn(E_Protein_maturation,...
                                          gene_ID,...
                                          sub_peptide,prod_peptide)
%peptide + h2o -> peptide-without-M + met               
RxnID = strcat('R_Met_aminopeptidase_',gene_ID);
Compo = {sub_peptide;'M_h2o_c';prod_peptide;'M_met__L_c'};
Coeff = [-1;-1;1;1];
Rev = 0;
Subproc = 'Protein maturation';
Catalyst = 'Methionine_Aminopeptidase_Enzyme_c';

E_Protein_maturation = MatrixOneReactionAdding(E_Protein_maturation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = NTerminusPrediction(gene,ID_new,ID_old,NT_IDnew,NT_lik,...
                                   NT_pos_processed)
if ismember(gene,NT_IDnew)
    index_NT = find(strcmp(gene,NT_IDnew));
    lik = NT_lik(index_NT);
    if lik < 80
    %pos of the case with likelyhood < 80% will be assumed to be 1
        pos = 1;
    else
        pos = NT_pos_processed(index_NT);
    end
else
    if ~ismember(gene,ID_old)
        pos = 1; %assumed to be 1
    else
        index_ID = find(strcmp(gene,ID_old));
        gene_newID = ID_new{index_ID};
        index_NT = find(strcmp(gene_newID,NT_IDnew));
        lik = NT_lik(index_NT);
        if lik < 80
        %pos of the case with likelyhood < 80% will be assumed to be 1
            pos = 1;
        else
            pos = NT_pos_processed(index_NT);
        end    
    end
end
end
