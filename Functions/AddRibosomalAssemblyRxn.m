%% AddRibosomalAssemblyRxn
%   This function will generate ribosomal assembly reactions representing
%   transcription elongation and termination.
function [E_Ribosomal_assembly,...
          E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = AddRibosomalAssemblyRxn(E_Ribosomal_assembly,...
                                                           Functional_protein,...
                                                           E_Enzyme_formation,...
                                                           E_Enzyme_dilution,...
                                                           E_Protein_degradation)

% Formation of ribosome_50S.
%formation of the catalyst ('llmg_0519', trigger factor)
RxnID = 'R_Ribosome_50S_Enzyme';
Substrate = FindFunctionalProtein('llmg_0519',Functional_protein);
Product = 'Ribosome_50S_Enzyme_c';
Compo = {Substrate;Product};
Coeff = [-1;1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 
%add protein degradation reactions
RxnID = 'R_Ribosome_50S_Enzyme_degradation';
Sub = strcat(Substrate(1:end-2),'_degradation_c');
Compo = {Product;Sub};
Coeff = [-1;1];
Rev = 0;
Subproc = 'Protein degradation';
Catalyst = '';
E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 
%add dilution reaction
RxnID = 'R_dilution_Ribosome_50S_Enzyme';
Compo = {Product};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 
%formation of ribosome_50S
RxnID = 'R_ribosome_50S';
Substrate = {'ribosome_50S_protein_c';...
             'rRNA_5S_unmodified_c';...
             'rRNA_23S_c'};
n_sub = length(Substrate);
Product = 'ribosome_50S_c';
Compo = cat(1,Substrate,Product);
Coeff = [-1*ones(n_sub,1);1];
Rev = 0;
Subproc = 'Ribosomal assembly';
Catalyst = 'Ribosome_50S_Enzyme_c';

E_Ribosomal_assembly = MatrixOneReactionAdding(E_Ribosomal_assembly,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 

% Formation of ribosome_30S.
%formation of the catalyst ('llmg_0371', Era; 'llmg_1791', RbfA;
%                           'llmg_0936', RimM)
RxnID = 'R_Ribosome_30S_Enzyme';
Era = FindFunctionalProtein('llmg_0371',Functional_protein);
RbfA = FindFunctionalProtein('llmg_1791',Functional_protein);
RimM = FindFunctionalProtein('llmg_0936',Functional_protein);
Substrate = {Era;RbfA;RimM};
n_sub = length(Substrate);
Product = 'Ribosome_30S_Enzyme_c';
Compo = cat(1,Substrate,Product);
Coeff = [-1*ones(n_sub,1);1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

%add protein degradation reactions
RxnID = 'R_Ribosome_30S_Enzyme_degradation';
Sub = cellfun(@(x) strcat(x(1:end-2),'_degradation_c'),...
                      Substrate,'UniformOutput',false);
Compo = cat(1,Product,Sub);
Coeff = [-1;-1*Coeff(1:end-1)];
Rev = 0;
Subproc = 'Protein degradation';
Catalyst = '';
E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 
%add dilution reaction
RxnID = 'R_dilution_Ribosome_30S_Enzyme';
Compo = {Product};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

%formation of ribosome_30S
RxnID = 'R_ribosome_30S';
Substrate = {'ribosome_30S_protein_c';'rRNA_16S_c';'M_gtp_c'};
n_sub = length(Substrate);
Product = {'ribosome_30S_c';'M_gdp_c';'M_pi_c';'M_h_c'};
n_prod = length(Product);
Compo = cat(1,Substrate,Product);
Coeff = [-1*ones(n_sub,1);1*ones(n_prod,1)];
Rev = 0;
Subproc = 'Ribosomal assembly';
Catalyst = 'Ribosome_30S_Enzyme_c';

E_Ribosomal_assembly = MatrixOneReactionAdding(E_Ribosomal_assembly,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

                                
% Formation of ribosome_70S.
RxnID = 'R_ribosome_70S';
Substrate = {'ribosome_30S_c';'ribosome_50S_c'};
Product = {'ribosome_70S_c'};
Compo = cat(1,Substrate,Product);
Coeff = [-1;-1;1];
Rev = 0;
Subproc = 'Ribosomal assembly';
Catalyst = '';

E_Ribosomal_assembly = MatrixOneReactionAdding(E_Ribosomal_assembly,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

% Degradation of ribosome_70S.    
RxnID = 'R_ribosome_70S_degradation';
Substrate = {'ribosome_70S_c'};
Product = {'ribosome_30S_protein_degraded_c';'ribosome_50S_protein_degraded_c';
           'rRNA_23S_c';'rRNA_16S_c';'rRNA_5S_unmodified_c'};
Compo = cat(1,Substrate,Product);
Coeff = [-1;1;1;1;1;1];
Rev = 0;
Subproc = 'Protein degradation';
Catalyst = '';
E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

%Dilution of ribosome_70S.
RxnID = 'R_dilution_ribosome_70S';
Compo = {'ribosome_70S_c'};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
