%% AddtRNAChargingRxn
%   This function will generate tRNA charging reactions.
%   The reaction formulation for this subprocess is based on KEGG pathway
%   database, i.e., Aminoacyl-tRNA biosynthesis.
%   http://www.kegg.jp/kegg-bin/highlight_pathway?scale=1.0&map=llm00970&keyword=
function [E_tRNA_charging,...
          E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = AddtRNAChargingRxn(E_Generic_RNA_renaming,...
                                                      Functional_protein,...
                                                      E_Enzyme_formation,...
                                                      E_Enzyme_dilution,...
                                                      E_Protein_degradation)

% Generate a list of generic tRNA.
index = E_Generic_RNA_renaming.CoeffList == 1;
Generic_RNA = E_Generic_RNA_renaming.CompoList(index);
index_tRNA = ~cellfun('isempty',(cellfun(@(x) strfind(x,'tRNA'),...
                  Generic_RNA,'UniformOutput',false)));
Generic_tRNA = unique(Generic_RNA(index_tRNA));

% Generate matrix.
E_tRNA_charging = MatrixGeneration;

%AA metabolite ID and abbreviation
cd Data;
AA_ID = importdata('AA_ID.txt');
cd ../;
AA_ID = split(AA_ID);
AA = struct();
AA.abbr_3 = AA_ID(:,1);
AA.metid = AA_ID(:,2);
AA.abbr_1 = AA_ID(:,3);

% charge Glu
AAstr = 'Glu';
Cata = 'llmg_2332';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Gln
%L. lactis does not have glutamine-tRNA synthetase, so it firstly produces
%tRNA(Gln)-Glu, and then tRNA(Gln)-Gln:
%1> tRNA(Gln) + Glu + ATP + H2O = tRNA(Gln)-Glu + AMP + PPi + H
%2> tRNA(Gln)-Glu + ATP + Gln = tRNA(Gln)-Gln + Glu + ADP + Pi
%The two steps is lumped in the model:
% tRNA(Gln) + Gln + 2 ATP + H2O = tRNA(Gln)-Gln + ADP + Pi + AMP + PPi + H
%Then the catalyst is assumed as a complex consists of glutamyl-tRNA
%synthetase and aspartyl-tRNA(Asn)/glutamyl-tRNA(Gln) amidotransferase.
AAstr = 'Gln';
Cata = 'llmg_2332 and llmg_0174 and llmg_0175 and llmg_0176';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingGlnRxn(E_tRNA_charging,...
                                            Functional_protein,...
                                            AAstr,AA,Generic_tRNA,Cata,...
                                            E_Enzyme_formation,...
                                            E_Enzyme_dilution,...
                                            E_Protein_degradation);

% charge Ala
AAstr = 'Ala';
Cata = 'llmg_1906';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);
                                
% charge Asp
AAstr = 'Asp';
Cata = 'llmg_2215';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Asn
AAstr = 'Asn';
Cata = 'llmg_2017';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Gly
AAstr = 'Gly';
Cata = 'llmg_1477 and llmg_1478';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Thr
AAstr = 'Thr';
Cata = 'llmg_2169';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Ser
AAstr = 'Ser';
Cata = 'llmg_0722';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Cys
AAstr = 'Cys';
Cata = 'llmg_2040';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Val
AAstr = 'Val';
Cata = 'llmg_2455';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Leu
AAstr = 'Leu';
Cata = 'llmg_1741';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Ile
AAstr = 'Ile';
Cata = 'llmg_2053';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Lys
AAstr = 'Lys';
Cata = 'llmg_0389';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Arg
AAstr = 'Arg';
Cata = 'llmg_2314';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Pro
AAstr = 'Pro';
Cata = 'llmg_2412';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);
% charge His
AAstr = 'His';
Cata = 'llmg_2217';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Phe
AAstr = 'Phe';
Cata = 'llmg_2195 and llmg_2196';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Tyr
AAstr = 'Tyr';
Cata = 'llmg_0401';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Trp
AAstr = 'Trp';
Cata = 'llmg_0079';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge Met
AAstr = 'Met';
Cata = 'llmg_1764';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                         Functional_protein,...
                                         AAstr,AA,Generic_tRNA,Cata,...
                                         E_Enzyme_formation,...
                                         E_Enzyme_dilution,...
                                         E_Protein_degradation);

% charge fMet
%1> tRNA(fMet) + Met + ATP + H2O = tRNA(fMet)-Met + AMP + PPi + H
%2> tRNA(fMet)-Met + 10fthf + H2O = tRNA(fMet)-fMet + thf + H
%The two steps is lumped in the model:
% tRNA(fMet) + Met + ATP + 2 H2O + 10fthf = tRNA(fMet)-fMet + thf + 2 H + AMP + PPi
%Then the catalyst is assumed as a complex consists of methionyl-tRNA
%synthetase (llmg_1764) and methionyl-tRNA formyltransferase (llmg_2147).
AAstr = 'fMet';
Cata = 'llmg_1764 and llmg_2147';
[E_tRNA_charging,...
 E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = addChargingfMetRxn(E_tRNA_charging,...
                                             Functional_protein,...
                                             AAstr,AA,Generic_tRNA,Cata,...
                                             E_Enzyme_formation,...
                                             E_Enzyme_dilution,...
                                             E_Protein_degradation);

%Charging for the codons that cannot be read by specific tRNA in L. lactis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that the case exists that there are two or more tRNA can read the
% codon that do not have its specific tRNA. For this case, the model only
% choose the tRNA with the lowest cost for its production.

cd Data;
[~,raw,~] = xlsread('AA&codon&tRNA.xlsx','not exist');
cd ../;
list_tRNA = raw(2:end,3);
list_tRNA_used = raw(2:end,6);

for i = 1:length(list_tRNA)
    tRNA_notexist = list_tRNA{i};
    tRNA_used = list_tRNA_used{i};
    
    RxnID = strcat('R_',tRNA_notexist(1:end-2));
    RxnID = strrep(RxnID,'charged_in','charging');
    substrate = {tRNA_used};
    product = {tRNA_notexist};
    Compo = cat(1,substrate,product);
    Coeff = [-1;1];
    Rev = 0;
    Subproc = 'tRNA charging';
    Catalyst = '';

    E_tRNA_charging = MatrixOneReactionAdding(E_tRNA_charging,...
                                                  RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                 
function [E_tRNA_charging,...
          E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = addChargingRxn(E_tRNA_charging,...
                                                  Functional_protein,...
                                                  AAstr,AA,Generic_tRNA,Cata,...
                                                  E_Enzyme_formation,...
                                                  E_Enzyme_dilution,...
                                                  E_Protein_degradation)
                                    
%AA metabolite ID and abbreviation
[~, index_AA] = ismember(AAstr,AA.abbr_3);
AA_met = AA.metid{index_AA};
AA_abbr = AA.abbr_1{index_AA};

%look for tRNA
index_tRNA = ~cellfun('isempty',...
                          (cellfun(@(x) strfind(x,strcat('_',AAstr)),...
                           Generic_tRNA,'UniformOutput',false)));
tRNAlist = Generic_tRNA(index_tRNA);
n_tRNA = length(tRNAlist);

for i = 1:n_tRNA
    tRNA = tRNAlist{i};
    
    %codon
    codon = tRNA(end-4:end-2);
    
    %add charging reaction
    RxnID = strcat('R_',AA_abbr,'_charging_',tRNA(1:end-2));
    Substrate = {AA_met;'M_atp_c';'M_h2o_c'};
    Product = {strcat(AA_abbr,'_charged_in_',tRNA(1:end-2),'_c');
               'M_amp_c';'M_ppi_c';'M_h_c'};
    Compo = cat(1,Substrate,Product);
    Coeff = [-1;-1;-1;1;1;1;1];
    Rev = 0;
    Subproc = 'tRNA charging';
    Enzyme = strcat('tRNA_Synthetase_',AA_abbr,'_',codon,'_Enzyme_c');
    Catalyst = Enzyme;

    E_tRNA_charging = MatrixOneReactionAdding(E_tRNA_charging,...
                                         RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
    %add catalyst formation reaction
    if isempty(strfind(Cata,'or')) && isempty(strfind(Cata,'and'))%if is a single enzyme
        RxnID = strcat('R_',Enzyme(1:end-2));
        Substrate = FindFunctionalProtein(Cata,Functional_protein);
        Product = Enzyme;
        Compo = {Substrate;Product};
        Coeff = [-1;1];
        Rev = 0;
        Subproc = 'Enzyme formation';
        Catalyst = '';

        E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
    
        %add protein degradation reactions
        RxnID = strcat('R_',Enzyme(1:end-2),'_degradation');
        Sub = strcat(Substrate(1:end-2),'_degradation_c');
        Compo = {Product;Sub};
        Coeff = [-1;1];
        Rev = 0;
        Subproc = 'Protein degradation';
        Catalyst = '';
        E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                         
        %add dilution reactions
        RxnID = strcat('R_dilution_',Enzyme(1:end-2));
        Compo = {Product};
        Coeff = -1;
        Rev = 0;
        Subproc = 'Enzyme dilution';
        Catalyst = '';
        E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);   
    
    elseif ~isempty(strfind(Cata,'and'))%if is a complex 
        Cata = strrep(Cata,' and ',' ');
        Cata = split(Cata);
        n_cata = length(Cata);
            
        RxnID = strcat('R_',Enzyme(1:end-2));            
        Substrate = cellfun(@(x) FindFunctionalProtein(x,Functional_protein),...
                                Cata,'UniformOutput',false);
        Product = {Enzyme};
        Compo = cat(1,Substrate,Product);
        Coeff = [-1*ones(n_cata,1);1];
        Rev = 0;
        Subproc = 'Enzyme formation';
        Catalyst = '';

        E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                                 RxnID,Compo,Coeff,Rev,Subproc,Catalyst);       

        %add protein degradation reaction
        RxnID = strcat('R_',Enzyme(1:end-2),'_degradation');
        Sub_raw = Compo(1:end-1);
        Sub = cellfun(@(x) strcat(x(1:end-2),'_degradation_c'),...
                              Sub_raw,'UniformOutput',false);
        Compo = cat(1,Product,Sub);
        Coeff = [-1;-1*Coeff(1:end-1)];
        Rev = 0;
        Subproc = 'Protein degradation';
        Catalyst = '';
        E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

        %add dilution reaction
        RxnID = strcat('R_dilution_',Enzyme(1:end-2));
        Compo = Product;
        Coeff = -1;
        Rev = 0;
        Subproc = 'Enzyme dilution';
        Catalyst = '';
        E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
        
	elseif ~isempty(strfind(Cata,'or'))%if has isozymes
        Cata = strrep(Cata,' or ',' ');
        Cata = split(Cata);
        n_cata = length(Cata);
        for j = 1:n_cata
            RxnID = strcat('R_',Enzyme(1:end-2),'_',num2str(j));
            Substrate = Cata{j};
            Substrate = FindFunctionalProtein(Substrate,Functional_protein);
            Product = Enzyme;      
            Compo = {Substrate;Product};
            Coeff = [-1;1];
            Rev = 0;
            Subproc = 'Enzyme formation';
            Catalyst = '';

            E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);            

            %add protein degradation reactions
            RxnID = strcat('R_',Enzyme(1:end-2),'_',num2str(j),'_degradation');
            Sub = strcat(Substrate(1:end-2),'_degradation_c');
            Compo = {Product;Sub};
            Coeff = [-1;1];
            Rev = 0;
            Subproc = 'Protein degradation';
            Catalyst = '';
            E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                                 RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

            %add dilution reactions
            RxnID = strcat('R_dilution_',Enzyme(1:end-2),'_',num2str(j));
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
            
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                 
function [E_tRNA_charging,...
          E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = addChargingGlnRxn(E_tRNA_charging,...
                                                     Functional_protein,...
                                                     AAstr,AA,Generic_tRNA,Cata,...
                                                     E_Enzyme_formation,...
                                                     E_Enzyme_dilution,...
                                                     E_Protein_degradation)
%AA metabolite ID and abbreviation
[~, index_AA] = ismember(AAstr,AA.abbr_3);
AA_met = AA.metid{index_AA};
AA_abbr = AA.abbr_1{index_AA};

%look for tRNA
index_tRNA = ~cellfun('isempty',...
                      (cellfun(@(x) strfind(x,strcat('_',AAstr)),...
                       Generic_tRNA,'UniformOutput',false)));
tRNAlist = Generic_tRNA(index_tRNA);
n_tRNA = length(tRNAlist);
for i = 1:n_tRNA
    tRNA = tRNAlist{i};
    
    %codon
    codon = tRNA(end-4:end-2);
    
    %add charging reaction
    %tRNA(Gln) + Gln + 2 ATP + H2O = tRNA(Gln)-Gln + ADP + Pi + AMP + PPi + H
    RxnID = strcat('R_',AA_abbr,'_charging_',tRNA(1:end-2));
    Substrate = {AA_met;'M_atp_c';'M_h2o_c'};
    Product = {strcat(AA_abbr,'_charged_in_',tRNA(1:end-2),'_c');
               'M_amp_c';'M_ppi_c';'M_h_c';'M_adp_c';'M_pi_c'};
    Compo = cat(1,Substrate,Product);
    Coeff = [-1;-2;-1;1;1;1;1;1;1];
    Rev = 0;
    Subproc = 'tRNA charging';
    Enzyme = strcat('tRNA_Synthetase_Complex_',AA_abbr,'_',codon,'_Enzyme_c');
    Catalyst = Enzyme;

    E_tRNA_charging = MatrixOneReactionAdding(E_tRNA_charging,...
                                         RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
    
    %add catalyst formation reaction
    Cata = strrep(Cata,' and ',' ');
    Cata = split(Cata);
    n_cata = length(Cata);
            
    RxnID = strcat('R_',Enzyme(1:end-2));            
    Substrate = cellfun(@(x) FindFunctionalProtein(x,Functional_protein),...
                                Cata,'UniformOutput',false);
    Product = {Enzyme};
    Compo = cat(1,Substrate,Product);
    Coeff = [-1*ones(n_cata,1);1];
    Rev = 0;
    Subproc = 'Enzyme formation';
    Catalyst = '';

    E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                              RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

    %add protein degradation reaction
    RxnID = strcat('R_',Enzyme(1:end-2),'_degradation');
    Sub_raw = Compo(1:end-1);
    Sub = cellfun(@(x) strcat(x(1:end-2),'_degradation_c'),...
                  Sub_raw,'UniformOutput',false);
    Compo = cat(1,Product,Sub);
    Coeff = [-1;-1*Coeff(1:end-1)];
    Rev = 0;
    Subproc = 'Protein degradation';
    Catalyst = '';
    E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

    %add dilution reaction
    RxnID = strcat('R_dilution_',Enzyme(1:end-2));
    Compo = Product;
    Coeff = -1;
    Rev = 0;
    Subproc = 'Enzyme dilution';
    Catalyst = '';
    E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                          
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                 
function [E_tRNA_charging,...
          E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = addChargingfMetRxn(E_tRNA_charging,...
                                                      Functional_protein,...
                                                      AAstr,AA,Generic_tRNA,Cata,...
                                                      E_Enzyme_formation,...
                                                      E_Enzyme_dilution,...
                                                      E_Protein_degradation)
%AA metabolite ID and abbreviation
[~, index_AA] = ismember('Met',AA.abbr_3);
AA_met = AA.metid{index_AA};
AA_abbr = AA.abbr_1{index_AA};

%look for tRNA
index_tRNA = ~cellfun('isempty',...
                          (cellfun(@(x) strfind(x,strcat('_',AAstr)),...
                           Generic_tRNA,'UniformOutput',false)));
tRNAlist = Generic_tRNA(index_tRNA);
n_tRNA = length(tRNAlist);
for i = 1:n_tRNA
    tRNA = tRNAlist{i};
    
    %codon
    codon = tRNA(end-4:end-2);
    
    %add charging reaction
    %tRNA(fMet) + Met + ATP + 2 H2O + 10fthf = tRNA(fMet)-fMet + thf + 2 H + AMP + PPi
    RxnID = strcat('R_f',AA_abbr,'_charging_',tRNA(1:end-2));
    Substrate = {AA_met;'M_atp_c';'M_h2o_c';'M_10fthf_c'};
    Product = {strcat('f',AA_abbr,'_charged_in_',tRNA(1:end-2),'_c');
               'M_amp_c';'M_ppi_c';'M_thf_c';'M_h_c'};
    Compo = cat(1,Substrate,Product);
    Coeff = [-1;-1;-2;-1;1;1;1;1;2];
    Rev = 0;
    Subproc = 'tRNA charging';
    Enzyme = strcat('tRNA_Synthetase_Complex_f',AA_abbr,'_',codon,'_Enzyme_c');
    Catalyst = Enzyme;

    E_tRNA_charging = MatrixOneReactionAdding(E_tRNA_charging,...
                                         RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
    
    %add catalyst formation reaction                                   
    Cata = strrep(Cata,' and ',' ');
    Cata = split(Cata);
    n_cata = length(Cata);
            
    RxnID = strcat('R_',Enzyme(1:end-2));            
    Substrate = cellfun(@(x) FindFunctionalProtein(x,Functional_protein),...
                                Cata,'UniformOutput',false);
    Product = {Enzyme};
    Compo = cat(1,Substrate,Product);
    Coeff = [-1*ones(n_cata,1);1];
    Rev = 0;
    Subproc = 'Enzyme formation';
    Catalyst = '';

    E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                              RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

    %add protein degradation reaction
    RxnID = strcat('R_',Enzyme(1:end-2),'_degradation');
    Sub_raw = Compo(1:end-1);
    Sub = cellfun(@(x) strcat(x(1:end-2),'_degradation_c'),...
                  Sub_raw,'UniformOutput',false);
    Compo = cat(1,Product,Sub);
    Coeff = [-1;-1*Coeff(1:end-1)];
    Rev = 0;
    Subproc = 'Protein degradation';
    Catalyst = '';
    E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

    %add dilution reaction
    RxnID = strcat('R_dilution_',Enzyme(1:end-2));
    Compo = Product;
    Coeff = -1;
    Rev = 0;
    Subproc = 'Enzyme dilution';
    Catalyst = '';
    E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

end
end
