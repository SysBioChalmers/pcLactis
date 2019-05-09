%% AddRibosomeProteinRxn
%   This function will generate ribosome protein formation reactions.
function [E_Ribosomal_assembly,...
          E_Protein_degradation]= AddRibosomeProteinRxn(E_Ribosomal_assembly,...
                                                        Functional_protein,...
                                                        E_Protein_degradation)

% Import ribosomal protein information.
cd Data;
[~, ~, RibLSUinfo] = xlsread('Ribosomal_protein.xlsx','LSU');
[~, ~, RibSSUinfo] = xlsread('Ribosomal_protein.xlsx','SSU');
cd ../;

% Formation and degradation of ribosome_50S_protein.
[E_Ribosomal_assembly,...
 E_Protein_degradation] = AddRibProtRxn(E_Ribosomal_assembly,...
                                        RibLSUinfo,...
                                        Functional_protein,...
                                        E_Protein_degradation);
% Formation and degradation of ribosome_30S_protein.                                         
[E_Ribosomal_assembly,...
 E_Protein_degradation] = AddRibProtRxn(E_Ribosomal_assembly,...
                                        RibSSUinfo,...
                                        Functional_protein,...
                                        E_Protein_degradation);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E_Ribosomal_assembly,...
          E_Protein_degradation] = AddRibProtRxn(E_Ribosomal_assembly,...
                                                 RibSUinfo,...
                                                 Functional_protein,...
                                                 E_Protein_degradation)
% Generate ribosomal protein formation reactions.
geneList = RibSUinfo(2:end,1);
subunitList = RibSUinfo(2:end,2);
coefficient = cell2mat(RibSUinfo(2:end,3));

% Formation.
if strcmp(subunitList{1}(1),'L')
    RxnID = 'R_ribosome_50S_protein';
    Product = 'ribosome_50S_protein_c';
elseif strcmp(subunitList{1}(1),'S')
    RxnID = 'R_ribosome_30S_protein';
    Product = 'ribosome_30S_protein_c';
end

Substrate = cellfun(@(x) FindFunctionalProtein(x,Functional_protein),...
                      geneList,'UniformOutput',false);
Compo = cat(1,Substrate,{Product});
Coeff = [-1*coefficient;1];
Rev = 0;
Subproc = 'Ribosomal assembly';
Catalyst = '';

E_Ribosomal_assembly = MatrixOneReactionAdding(E_Ribosomal_assembly,...
                                 RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

% Degradation into subunits.
RxnID = strcat(RxnID,'_degraded_degradation');
Sub = strcat(Product(1:end-2),'_degraded_c');
Prod = cellfun(@(x) strcat(x(1:end-2),'_degradation_c'),...
                      Substrate,'UniformOutput',false);
Compo = cat(1,{Sub},Prod);
Coeff = [-1;coefficient];
Rev = 0;
Subproc = 'Protein degradation';
Catalyst = '';
E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);



end