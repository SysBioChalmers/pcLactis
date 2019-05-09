%% AddTranslationRxn
%   This function will generate translation reactions.
function [E_Translation,...
          Translated_gene,...
          E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = AddTranslationRxn(TU_degraded,...
                                                     TU,Sub_TU_cleaved,Gene_list,...
                                                     Seq_gene_NA_total,...
                                                     Seq_CDS_AA_total,...
                                                     E_tRNA_charging,...
                                                     Functional_protein,...
                                                     E_Enzyme_formation,...
                                                     E_Enzyme_dilution,...
                                                     E_Protein_degradation)

% Sum of TUs and Sub-TUs.
TU_translation_pool = struct();
TU_translation_pool.ID = TU_degraded;
TU_translation_pool.sequence = {};
for i = 1:length(TU_translation_pool.ID)
    TU_name = TU_translation_pool.ID{i};
    if ismember(TU_name,TU.TUID)
        [~, index_TU] = ismember(TU_name,TU.TUID);
        TU_seq = TU.sequence{index_TU};
        TU_translation_pool.sequence = cat(1,TU_translation_pool.sequence,{TU_seq});
    elseif ismember(TU_name,Sub_TU_cleaved.ID)
        [~, index_TU] = ismember(TU_name,Sub_TU_cleaved.ID);
        TU_seq = Sub_TU_cleaved.sequence{index_TU};
        TU_translation_pool.sequence = cat(1,TU_translation_pool.sequence,{TU_seq});
    end
end

% List of proteins can be translated in the model, which is the same as the
% list of Gene_list.total_protein_gene.
Translated_gene = struct();
Translated_gene.ID = Gene_list.total_protein_gene;
Translated_gene.NAseq = {};
Translated_gene.AAseq = {};
for i = 1:length(Translated_gene.ID)
    gene_name = Translated_gene.ID{i};
    gene_NAseq = GetSequence(gene_name,Seq_gene_NA_total.gene,Seq_gene_NA_total.sequence);
    gene_AAseq = GetSequence(gene_name,Seq_CDS_AA_total.gene,Seq_CDS_AA_total.sequence);
    Translated_gene.NAseq = cat(1,Translated_gene.NAseq,{gene_NAseq});
    Translated_gene.AAseq = cat(1,Translated_gene.AAseq,{gene_AAseq});
end

% Generate a list of charged generic tRNA.
index = E_tRNA_charging.CoeffList == 1;
Product = E_tRNA_charging.CompoList(index);
index_charged_tRNA = ~cellfun('isempty',(cellfun(@(x) strfind(x,'tRNA'),...
                              Product,'UniformOutput',false)));
Charged_tRNA = unique(Product(index_charged_tRNA));

% Generate matrix.
E_Translation = MatrixGeneration;

% Translation initiation factor complex.
[E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = AddIFComplexRxn(E_Enzyme_formation,...
                                          E_Enzyme_dilution,...
                                          E_Protein_degradation,...
                                          Functional_protein);

% EF-Tu-EF-Ts-complex.
[E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = AddTuTsComplexRxn(E_Enzyme_formation,...
                                            E_Enzyme_dilution,...
                                            E_Protein_degradation,...
                                            Functional_protein);

% EF-G.
[E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = AddEFGRxn(E_Enzyme_formation,...
                                    E_Enzyme_dilution,...
                                    E_Protein_degradation,...
                                    Functional_protein);

% Translation termination complex.
[E_Enzyme_formation,...
 E_Enzyme_dilution,...
 E_Protein_degradation] = AddRFRrfRxn(E_Enzyme_formation,...
                                      E_Enzyme_dilution,...
                                      E_Protein_degradation,...
                                      Functional_protein);

% Add tRNA activation step:
%A-tRNA-Ala + gtp + h2o -> A-tRNA-Ala-activated + gdp + pi + h
%The model integrates EF-Tu and EF-Ts into a complex as the catalyst of the
%reaction.
index = cellfun('isempty',cellfun(@(x) strfind(x,'fMet'),...
                                       Charged_tRNA,...
                                       'UniformOutput',false));
Charged_tRNA_for_elongation = Charged_tRNA(index);%remove tRNA-fMet
E_Translation = AddtRNAActivationStep(E_Translation,...
                                      Charged_tRNA_for_elongation);

% Add tRNA elongation step:
%A-tRNA-Ala-activated + gtp + h2o -> A-tRNA-Ala-elongated + gdp + pi + h
%The substrate of this reaction is used as a unit in translation elongation
%reaction. The catalyst of the reaction is TU-G.
E_Translation = AddtRNAElongationStep(E_Translation,...
                                      Charged_tRNA_for_elongation);

% Generate a structure containing relationship between codons and AAs.
DB = struct();
DB.subtRNA = {};
DB.codon = {};
DB.AA = {};
for i = 1:length(Charged_tRNA_for_elongation)
    raw = Charged_tRNA_for_elongation{i};
    DB.subtRNA = cat(1,DB.subtRNA,strcat(raw(1:end-2),'_elongated_c'));
    DB.codon = cat(1,DB.codon,raw(end-4:end-2));
    DB.AA = cat(1,DB.AA,raw(1));
end
DB.codon = cellfun(@(x) strrep(x,'U','T'),...%replace U with T in codon
                   DB.codon,'UniformOutput',false);
                                  
% Add translation reactions.
for i = 1:length(Translated_gene.ID)
    gene_ID = Translated_gene.ID{i};
    gene_NA = Translated_gene.NAseq{i};
    gene_AA = Translated_gene.AAseq{i};
    index_TU = cellfun(@(x) strfind(x,gene_NA),...
                       TU_translation_pool.sequence,...
                       'UniformOutput',false);
    index_TU = cellfun('isempty',index_TU) == 0;
    TU_list = TU_translation_pool.ID(index_TU);
    
    %add translation initiation reactions
    Prod_initiated = {};
    for j = 1:length(TU_list)%a gene could be from more than one TU
        TU_ID = TU_list{j};
        [E_Translation,...
         prod_init]= AddTranslationInitiationRxn(E_Translation,...
                                                 gene_ID,TU_ID);
        Prod_initiated = cat(1,Prod_initiated,prod_init);
    end
    
    %add translation elongation reactions
    Prod_elongated = {};
    for j = 1:length(Prod_initiated)%a gene could have more than one initiated product
        sub_ID = Prod_initiated{j};
        [E_Translation, prod_elong] = AddTranslationElongationRxn(E_Translation,...
                                                                  gene_ID,sub_ID,...
                                                                  DB,...
                                                                  gene_NA,gene_AA);    
        Prod_elongated = cat(1,Prod_elongated,prod_elong);
    end
    
     %add translation termination reactions
    for j = 1:length(Prod_elongated)%a gene could have more than one elongated product
        sub_ID = Prod_elongated{j};
        E_Translation = AddTranslationTerminationRxn(E_Translation,...
                                                     gene_ID,sub_ID,...
                                                     gene_NA);
    end
end

% Import to Excel file.
%model = Matrix2CBModel(E_Translation);
%writeCbModel(model,'xls','Translation.xls');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E_Translation,...
          prod_init]= AddTranslationInitiationRxn(E_Translation,...
                                                  gene_ID,TU_ID)
%translation initiation reaction
%fM-tRNA-fMet + h2o + gtp -> TU-ini + pi + h + gdp                                                 
RxnID = strcat('R_translation_initiation_',TU_ID,'_',gene_ID);
prod_init = strcat(TU_ID,'_',gene_ID,'_initiated_c');
Compo = {'fM_charged_in_tRNA_fMet_AUG_c';'M_h2o_c';'M_gtp_c';
         prod_init;'M_pi_c';'M_h_c';'M_gdp_c'};
Coeff = [-1;-1;-1;1;1;1;1];
Rev = 0;
Subproc = 'Translation';
Catalyst = 'Translation_Initiation_Complex_Enzyme_c';

E_Translation = MatrixOneReactionAdding(E_Translation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = AddIFComplexRxn(E_Enzyme_formation,...
                                                   E_Enzyme_dilution,...
                                                   E_Protein_degradation,...
                                                   Functional_protein)
%IF1 + IF2 + IF3 -> IF-complex
RxnID = 'R_Translation_Initiation_Complex_Enzyme';
IF1 = FindFunctionalProtein('llmg_2358',Functional_protein);
IF2 = FindFunctionalProtein('llmg_1792',Functional_protein);
IF3 = FindFunctionalProtein('llmg_2031',Functional_protein);
Product = 'Translation_Initiation_Complex_Enzyme_c';
Compo = {IF1;IF2;IF3;Product};
Coeff = [-1;-1;-1;1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 
%add protein degradation reactions
RxnID = 'R_Translation_Initiation_Complex_Enzyme_degradation';
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
RxnID = 'R_dilution_Translation_Initiation_Complex_Enzyme';
Compo = {Product};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_Translation = AddtRNAActivationStep(E_Translation,...
                                               Charged_tRNA_for_elongation)
%A-tRNA-Ala + gtp + h2o -> A-tRNA-Ala-activated + gdp + pi + h
for i = 1:length(Charged_tRNA_for_elongation)
    tRNA_activate = Charged_tRNA_for_elongation{i};
    tRNA_activate = tRNA_activate(1:end-2);
    %tRNA activation step
    RxnID = strcat('R_Activate_',tRNA_activate);
    Compo = {strcat(tRNA_activate,'_c');'M_gtp_c';'M_h2o_c';
             strcat(tRNA_activate,'_activated_c');'M_pi_c';'M_h_c';'M_gdp_c'};
    Coeff = [-1;-1;-1;1;1;1;1];
    Rev = 0;
    Subproc = 'Translation';
    Catalyst = 'EF_Tu_EF_Ts_Complex_Enzyme_c';

    E_Translation = MatrixOneReactionAdding(E_Translation,...
                                            RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = AddTuTsComplexRxn(E_Enzyme_formation,...
                                                     E_Enzyme_dilution,...
                                                     E_Protein_degradation,...
                                                     Functional_protein)
%EF-Tu + EF-Ts -> EF-Tu-EF-Ts-complex
RxnID = 'R_EF_Tu_EF_Ts_Complex_Enzyme';
EF_Tu = FindFunctionalProtein('llmg_2050',Functional_protein);
EF_Ts = FindFunctionalProtein('llmg_2429',Functional_protein);
Product = 'EF_Tu_EF_Ts_Complex_Enzyme_c';
Compo = {EF_Tu;EF_Ts;Product};
Coeff = [-1;-1;1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

%add protein degradation reactions
RxnID = 'R_EF_Tu_EF_Ts_Complex_Enzyme_degradation';
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
RxnID = 'R_dilution_EF_Tu_EF_Ts_Complex_Enzyme';
Compo = {Product};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_Translation = AddtRNAElongationStep(E_Translation,...
                                               Charged_tRNA_for_elongation)
%A-tRNA-Ala-activated + gtp + h2o -> A-tRNA-Ala-elongated + gdp + pi + h
for i = 1:length(Charged_tRNA_for_elongation)
    tRNA_elongate = Charged_tRNA_for_elongation{i};
    tRNA_elongate = tRNA_elongate(1:end-2);
    %tRNA elongation step
    RxnID = strcat('R_Elongate_',tRNA_elongate);
    Compo = {strcat(tRNA_elongate,'_activated_c');'M_gtp_c';'M_h2o_c';
             strcat(tRNA_elongate,'_elongated_c');'M_pi_c';'M_h_c';'M_gdp_c'};
    Coeff = [-1;-1;-1;1;1;1;1];
    Rev = 0;
    Subproc = 'Translation';
    Catalyst = 'EF_G_Enzyme_c';

    E_Translation = MatrixOneReactionAdding(E_Translation,...
                                            RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = AddEFGRxn(E_Enzyme_formation,...
                                             E_Enzyme_dilution,...
                                             E_Protein_degradation,...
                                             Functional_protein)
%EF-G -> EF-G
RxnID = 'R_EF_G_Enzyme';
EF_G = FindFunctionalProtein('llmg_2556',Functional_protein);
Product = 'EF_G_Enzyme_c';
Compo = {EF_G;Product};
Coeff = [-1;1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 
%add protein degradation reactions
RxnID = 'R_EF_G_Enzyme_degradation';
EF_G = strcat(EF_G(1:end-2),'_degradation_c');
Compo = {Product;EF_G};
Coeff = [-1;1];
Rev = 0;
Subproc = 'Protein degradation';
Catalyst = '';
E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                         
%add dilution reactions
RxnID = strcat('R_dilution_EF_G_Enzyme');
Compo = {Product};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E_Translation, prod_elong] = AddTranslationElongationRxn(E_Translation,...
                                                                  gene_ID,sub_ID,...
                                                                  DB,...
                                                                  gene_NA,gene_AA)
%translation elongation reaction
%TU-ini + A-tRNA-Ala-elongated + ... + h2o -> TU-elongated
RxnID = strcat('R_translation_elongation_',sub_ID(1:end-12));
prod_elong = strrep(sub_ID,'initiated','elongated');

%count numbers of each tRNA used
tRNA = struct();
tRNA.sub = DB.subtRNA;
tRNA.num =zeros(length(tRNA.sub),1);

%generate a structure to compare translated AA seq by the code with AA seq
%given by NCBI.
AA = struct();
AA.ID = unique(DB.AA);
AA.seq_trans = zeros(length(AA.ID),1);
AA.seq_NCBI = zeros(length(AA.ID),1);

%calculate the number of tRNA used
if ismember(gene_NA(end-2:end),{'TAG','TAA','TGA'})
    gene_NA = gene_NA(4:end-3);
else%some genes do not end with stop codon
    gene_NA = gene_NA(4:end);
end

for n = 1:length(gene_NA)/3
    triplet = gene_NA(3*n-2:3*n);
    [~, index] = ismember(triplet,DB.codon);
    tRNA.num(index) = tRNA.num(index)+1;
end

%check if the translated AA seq is in line with that from NCBI
for m = 1:length(tRNA.num)
    if tRNA.num(m) > 0
        AAinfo = tRNA.sub{m};
        AAabbr = AAinfo(1);
        [~, index_AA] = ismember(AAabbr,AA.ID);
        AA.seq_trans(index_AA) = AA.seq_trans(index_AA)+tRNA.num(m);
    end
end
gene_AA_check = gene_AA(2:end);
for m = 1:length(AA.ID)
	AA.seq_NCBI(m) = count(gene_AA_check,AA.ID{m});
end
if ~all(AA.seq_trans == AA.seq_NCBI)
	error(['The AA sequence of ',gene_ID,...
          ' translated is different from that from NCBI.']);
end

index_coef = find(tRNA.num);
tRNA_sub = tRNA.sub(index_coef);
coef_tRNA_sub = -1*tRNA.num(index_coef);

coef_h2o = length(gene_NA)/3;%equals to the number of AAs added.

Compo = cat(1,{sub_ID},tRNA_sub,{'M_h2o_c'},{prod_elong});
Coeff = [-1;coef_tRNA_sub;coef_h2o;1];
Rev = 0;
Subproc = 'Translation';
Catalyst = 'ribosome_70S_c';

E_Translation = MatrixOneReactionAdding(E_Translation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = AddRFRrfRxn(E_Enzyme_formation,...
                                               E_Enzyme_dilution,...
                                               E_Protein_degradation,...
                                               Functional_protein)

%RF1 + RF3 + Rrf -> RFRrf_UAA_UAG
RxnID = 'R_RFRrf_UAA_UAG_Enzyme';
RF1 = FindFunctionalProtein('llmg_0557',Functional_protein);
RF3 = FindFunctionalProtein('llmg_0368',Functional_protein);
Rrf = FindFunctionalProtein('llmg_2284',Functional_protein);
Product = 'RFRrf_UAA_UAG_Enzyme_c';
Compo = {RF1;RF3;Rrf;Product};
Coeff = [-1;-1;-1;1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

%add protein degradation reactions
RxnID = 'R_RFRrf_UAA_UAG_Enzyme_degradation';
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
RxnID = 'R_dilution_RFRrf_UAA_UAG_Enzyme';
Compo = {Product};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 
%RF2 + RF3 + Rrf -> RFRrf_UGA
RxnID = 'R_RFRrf_UGA_Enzyme';
RF2 = FindFunctionalProtein('llmg_1547',Functional_protein);
RF3 = FindFunctionalProtein('llmg_0368',Functional_protein);
Rrf = FindFunctionalProtein('llmg_2284',Functional_protein);
Product = 'RFRrf_UGA_Enzyme_c';
Compo = {RF2;RF3;Rrf;Product};
Coeff = [-1;-1;-1;1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

%add protein degradation reactions
RxnID = 'R_RFRrf_UGA_Enzyme_degradation';
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
RxnID = 'R_dilution_RFRrf_UGA_Enzyme';
Compo = {Product};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                             RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_Translation = AddTranslationTerminationRxn(E_Translation,...
                                                      gene_ID,sub_ID,...
                                                      gene_NA)
%translation termination reaction
%TU-elongated + gtp + h2o -> peptide + gdp + pi + h
RxnID = strcat('R_translation_termination_',sub_ID(1:end-12));
prod_termin = strcat(gene_ID,'_nascent_c');
Compo = {sub_ID;'M_gtp_c';'M_h2o_c';
         prod_termin;'M_pi_c';'M_h_c';'M_gdp_c'};
Coeff = [-1;-1;-1;1;1;1;1];
Rev = 0;
Subproc = 'Translation';

if ismember(gene_NA(end-2:end),{'TAG'})%end with 'TAG' then catalysed by RFRrf_UAG
    Catalyst = 'RFRrf_UAA_UAG_Enzyme_c';
elseif ismember(gene_NA(end-2:end),{'TGA'})%end with 'TGA' then catalysed by RFRrf_UGA
    Catalyst = 'RFRrf_UGA_Enzyme_c';
else%end with 'TAA' or without stop codon then catalysed by RFRrf_UAA
    Catalyst = 'RFRrf_UAA_UAG_Enzyme_c';
end

E_Translation = MatrixOneReactionAdding(E_Translation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
end
