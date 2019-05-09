%% AddProteinDegradationRxn
%   This function will generate protein degradation reactions.
function [E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = AddProteinDegradationRxn(Functional_protein,...
                                                            Translated_gene,...
                                                            E_Enzyme_formation,...
                                                            E_Enzyme_dilution,...
                                                            E_Protein_degradation)
% Add formation of protein degradation complex.

ClpP = FindFunctionalProtein('llmg_0638',Functional_protein);
ClpB = FindFunctionalProtein('llmg_0986',Functional_protein);
ClpC = FindFunctionalProtein('llmg_0615',Functional_protein);
ClpE = FindFunctionalProtein('llmg_0528',Functional_protein);
Product = 'Clp_Protease_Complex_Enzyme_c';

RxnID = 'R_Clp_Protease_Complex_Enzyme';
Compo = {ClpP;ClpB;ClpC;ClpE;Product};
Coeff = [-1;-1;-1;-1;1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                 RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

% Add protein degradation reaction.    
RxnID = 'R_Clp_Protease_Complex_Enzyme_degradation';
Sub_raw = Compo(1:end-1);
Sub = cellfun(@(x) strcat(x(1:end-2),'_degradation_c'),...
                      Sub_raw,'UniformOutput',false);
Compo = cat(1,{Product},Sub);
Coeff = [-1;-1*Coeff(1:end-1)];
Rev = 0;
Subproc = 'Protein degradation';
Catalyst = '';
E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

% Add dilution reaction.
RxnID = 'R_dilution_Clp_Protease_Complex_Enzyme';
Compo = {Product};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 
% Add reactions for protein degradation.

%collect proteins to be degraded
idx = contains(E_Protein_degradation.CompoList,'_degradation_c');
protein_degraded_list = unique(E_Protein_degradation.CompoList(idx));

%protein degradation reactions
for i = 1:length(protein_degraded_list)
    protein_degraded = protein_degraded_list{i};%degraded protein ID
    
    protein_with_stoich = strrep(protein_degraded,'_degradation_c','');%protein stoichiometry
    
    if contains(protein_with_stoich,'_assumed')
        protein_with_stoich = strrep(protein_with_stoich,'_assumed','');
    end
    
    idx = max(strfind(protein_with_stoich,'_'));
    protein_id = protein_with_stoich(1:idx-1);%protein ID
    
    [~, idx] = ismember(protein_id,Translated_gene.ID);
    protein_mature_seq = Translated_gene.AAseq_mature{idx};%mature protein sequence

    sum = CountAA(protein_mature_seq);%count number of AAs
    
    RxnID = strcat('R_',protein_with_stoich,'_degradation');
    Substrate = {'M_atp_c';'M_h2o_c';protein_degraded};
    Product = cat(1,{'M_pi_c';'M_h_c';'M_adp_c'},sum.metid);
    Compo = cat(1,Substrate,Product);
    Rev = 0;
    Subproc = 'Protein degradation';
    Catalyst = 'Clp_Protease_Complex_Enzyme_c';
    
    if contains(protein_with_stoich,'Monomer')%the protein is a monomer
        %determine coefficient
        nAA = length(protein_mature_seq);
        natp = round(nAA*0.25);%the factor can be adjusted
        nadp = natp;
        npi = natp;
        nh2o = natp+nAA-1;
        nh = nh2o;
        
        Coeff_sub = -1*[natp;nh2o;1];
        Coeff_prod = [npi;nh;nadp;sum.num];
        Coeff = [Coeff_sub;Coeff_prod];
        
        E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
    else%is a oligomer
        %get stoichiometry
        z = strfind(protein_with_stoich,'mer');
        a = max(strfind(protein_with_stoich(1:z),'_'));
        x = str2double(protein_with_stoich(a+1:z-1));
        %determine coefficient
        nAA = length(protein_mature_seq);
        natp = x*round(nAA*0.25);
        nadp = natp;
        npi = natp;
        nh2o = natp+x*(nAA-1);
        nh = nh2o;
        
        Coeff_sub = -1*[natp;nh2o;1];
        Coeff_prod = [npi;nh;nadp;x*sum.num];
        Coeff = [Coeff_sub;Coeff_prod];
        
        E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
    end
end
