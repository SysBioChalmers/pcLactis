%% AddtRNAModificationRxn
function [E_tRNA_modification,...
          FUtRNA,...
          E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = AddtRNAModificationRxn(Functional_protein,...
                                                          E_Enzyme_formation,...
                                                          E_Enzyme_dilution,...
                                                          E_Protein_degradation)

% Generate matrix.
E_tRNA_modification = MatrixGeneration;

% Import tRNA modification data.
cd Data;
[~, ~, tRNAmod] = xlsread('RNA_modification.xlsx','tRNA');
[~, ~, catalystlist] = xlsread('RNA_modification.xlsx','Catalysts');
cd ../;

% Add enzyme reactions.
[E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = AddCatalystRxn(E_Enzyme_formation,...
                                         E_Enzyme_dilution,catalystlist,...
                                         Functional_protein,...
                                         E_Protein_degradation);

% Add tRNA modification reactions.
data = struct();
data.gene = tRNAmod(2:end,1);
data.prod = tRNAmod(2:end,2);
%change tRNA product to codon reader
data.prod = cellfun(@(x) ...
    strcat(x(1:length(x)-3),reverse(RNASeqComplement(x(length(x)-2:end)))),...
    data.prod,'UniformOutput',false);
data.cata = tRNAmod(2:end,5);
data.rxn = tRNAmod(2:end,6);

unique_gene = unique(data.gene);

FUtRNA = {};%collect totally modified tRNA

for i = 1:length(unique_gene)
    gene = unique_gene{i};
    n_data = find(strcmp(data.gene,gene));
    n = length(n_data);
    
    if n == 1%only one modification site
        for j = 1
            RxnID = strcat('R_tRNA_modification_',gene,'_',data.prod{n_data(j)});
            tRNAsub = strcat(gene,'_unmodified_c');
            tRNAprod = strcat(gene,'_',data.prod{n_data(j)},'_c');

            E_tRNA_modification = AddtRNAMdfRxn(E_tRNA_modification,...
                                                RxnID,tRNAsub,tRNAprod,...
                                                data,n_data,j);

            FUtRNA = cat(1,FUtRNA,tRNAprod);
        end
        
    elseif n == 2%two modification sites
        
        for j = 1
            RxnID = strcat('R_tRNA_modification_',gene,'_modified_',num2str(j));
            tRNAsub = strcat(gene,'_unmodified_c');
            tRNAprod = strcat(gene,'_modified_',num2str(j),'_c');
            
            E_tRNA_modification = AddtRNAMdfRxn(E_tRNA_modification,...
                                                RxnID,tRNAsub,tRNAprod,...
                                                data,n_data,j);

        end
        for j = 2
            RxnID = strcat('R_tRNA_modification_',gene,'_',data.prod{n_data(j)});
            tRNAsub = strcat(gene,'_modified_',num2str(j-1),'_c');
            tRNAprod = strcat(gene,'_',data.prod{n_data(j)},'_c');
            
            E_tRNA_modification = AddtRNAMdfRxn(E_tRNA_modification,...
                                                RxnID,tRNAsub,tRNAprod,...
                                                data,n_data,j);
        
            FUtRNA = cat(1,FUtRNA,tRNAprod);
        end

    elseif n > 2%more than two modification sites
        for j = 1
            RxnID = strcat('R_tRNA_modification_',gene,'_modified_',num2str(j));
            tRNAsub = strcat(gene,'_unmodified_c');
            tRNAprod = strcat(gene,'_modified_',num2str(j),'_c');
            
            E_tRNA_modification = AddtRNAMdfRxn(E_tRNA_modification,...
                                                RxnID,tRNAsub,tRNAprod,...
                                                data,n_data,j);
        end
        
        for j = 2:n-1
            RxnID = strcat('R_tRNA_modification_',gene,'_modified_',num2str(j));
            tRNAsub = strcat(gene,'_modified_',num2str(j-1),'_c');
            tRNAprod = strcat(gene,'_modified_',num2str(j),'_c');
            
            E_tRNA_modification = AddtRNAMdfRxn(E_tRNA_modification,...
                                                RxnID,tRNAsub,tRNAprod,...
                                                data,n_data,j);
        end
        
        for j = n
            RxnID = strcat('R_tRNA_modification_',gene,'_',data.prod{n_data(j)});
            tRNAsub = strcat(gene,'_modified_',num2str(j-1),'_c');
            tRNAprod = strcat(gene,'_',data.prod{n_data(j)},'_c');
            
            E_tRNA_modification = AddtRNAMdfRxn(E_tRNA_modification,...
                                                RxnID,tRNAsub,tRNAprod,...
                                                data,n_data,j);
                                           
            FUtRNA = cat(1,FUtRNA,tRNAprod);
        end        
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_tRNA_modification = AddtRNAMdfRxn(E_tRNA_modification,...
                                             RxnID,tRNAsub,tRNAprod,...
                                             data,n_data,j)

if all(isnan(data.rxn{n_data(j)}))%no other metabolites involved
	if all(isnan(data.cata{n_data(j)}))%no catalyst involved
        Catalyst = '';
        E_tRNA_modification = AddRxnForNoMet(E_tRNA_modification,...
                                             RxnID,tRNAsub,tRNAprod,...
                                             Catalyst);
	else%catalyst involved
        Cata = data.cata{n_data(j)};
        if contains(Cata,'and')
            Cata = strrep(Cata,' and ','_');
        end
        %add tRNA modification reaction
        Catalyst = strcat('RNAmod_',Cata,'_Enzyme_c');
        E_tRNA_modification = AddRxnForNoMet(E_tRNA_modification,...
                                             RxnID,tRNAsub,tRNAprod,...
                                             Catalyst);
	end
else%other metabolites involved
	MetInfo = extractBetween(data.rxn{n_data(j)},'(',')');
	Metname = MetInfo{1};
    Metname = strrep(Metname,', ',' ');
    Metname = transpose(strsplit(Metname));    
    Metcoeff = transpose(str2num(MetInfo{2}));
            
	if all(isnan(data.cata{n_data(j)}))%no catalyst involved
        Catalyst = '';
        E_tRNA_modification = AddRxnForMet(E_tRNA_modification,...
                                           RxnID,tRNAsub,tRNAprod,...
                                           Metname,Metcoeff,Catalyst);
	else%catalyst involved
        Cata = data.cata{n_data(j)};
        if contains(Cata,'and')
            Cata = strrep(Cata,' and ','_');
        end
        %add tRNA modification reaction
        Catalyst = strcat('RNAmod_',Cata,'_Enzyme_c');
        E_tRNA_modification = AddRxnForMet(E_tRNA_modification,...
                                           RxnID,tRNAsub,tRNAprod,...
                                           Metname,Metcoeff,Catalyst);
	end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_tRNA_modification = AddRxnForNoMet(E_tRNA_modification,...
                                              RxnID,tRNAsub,tRNAprod,...
                                              Catalyst)
Compo = {tRNAsub;tRNAprod};
Coeff = [-1;1];
Rev = 0;
Subproc = 'tRNA modification';

E_tRNA_modification = MatrixOneReactionAdding(E_tRNA_modification,...
                                      RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_tRNA_modification = AddRxnForMet(E_tRNA_modification,...
                                            RxnID,tRNAsub,tRNAprod,...
                                            Metname,Metcoeff,Catalyst)
Compo = cat(1,{tRNAsub;tRNAprod},Metname);
Coeff = [-1;1;Metcoeff];
Rev = 0;
Subproc = 'tRNA modification';

E_tRNA_modification = MatrixOneReactionAdding(E_tRNA_modification,...
                                      RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = AddCatalystRxn(E_Enzyme_formation,...
                                                  E_Enzyme_dilution,catalystlist,...
                                                  Functional_protein,...
                                                  E_Protein_degradation)
n = length(catalystlist);
for i = 1:n
    Cata = catalystlist{i};
    if ~contains(Cata,'and')%single enzyme

        Substrate = FindFunctionalProtein(Cata,Functional_protein);
        Product = strcat('RNAmod_',Cata,'_Enzyme_c');
        
        %add enzyme formation reactions
        RxnID = strcat('R_RNAmod_',Cata,'_Enzyme');
        Compo = {Substrate;Product};
        Coeff = [-1;1];
        Rev = 0;
        Subproc = 'Enzyme formation';
        Catalyst = '';
        E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                         RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                     
        %add protein degradation reactions
        RxnID = strcat('R_RNAmod_',Cata,'_Enzyme_degradation');
        Sub = strcat(Substrate(1:end-2),'_degradation_c');
        Compo = {Product;Sub};
        Coeff = [-1;1];
        Rev = 0;
        Subproc = 'Protein degradation';
        Catalyst = '';
        E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                         
        %add dilution reactions
        RxnID = strcat('R_dilution_RNAmod_',Cata,'_Enzyme');
        Compo = {Product};
        Coeff = -1;
        Rev = 0;
        Subproc = 'Enzyme dilution';
        Catalyst = '';
        E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                     
    else%a complex containing more than one subunit
        Cataname = strrep(Cata,' and ','_');
        Cata = strrep(Cata,' and ',' ');
        Cata = transpose(strsplit(Cata));
        Substrate = cellfun(@(x) FindFunctionalProtein(x,Functional_protein),...
                       Cata,'UniformOutput',false);
        nCata = length(Substrate);
        Product = strcat('RNAmod_',Cataname,'_Enzyme_c');
        
        %add enzyme formation reactions
        RxnID = strcat('R_RNAmod_',Cataname,'_Enzyme');
        Compo = cat(1,Substrate,Product);
        Coeff = [-1*ones(nCata,1);1];
        Rev = 0;
        Subproc = 'Enzyme formation';
        Catalyst = '';
        E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                         RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                     
        %add protein degradation reactions   
        RxnID = strcat('R_RNAmod_',Cataname,'_Enzyme_degradation');
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
                                     
        %add dilution reactions
        RxnID = strcat('R_dilution_RNAmod_',Cataname,'_Enzyme');
        Compo = {Product};
        Coeff = -1;
        Rev = 0;
        Subproc = 'Enzyme dilution';
        Catalyst = '';
        E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
    end
end
end
