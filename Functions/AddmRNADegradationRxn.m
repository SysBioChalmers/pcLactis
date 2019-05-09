%% AddmRNADegradationRxn
%   This function will generate mRNA degradation reactions for TUs and
%   Sub-TUs without stable RNA genes.
function [E_mRNA_degradation,...
          TU_degraded,...
          E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation,...
          E_RNA_dilution] = AddmRNADegradationRxn(Gene_list,TU_cleaved,...
                                                  Sub_TU_cleaved,...
                                                  TU,Functional_protein,...
                                                  E_Enzyme_formation,...
                                                  E_Enzyme_dilution,...
                                                  E_Protein_degradation,...
                                                  E_RNA_dilution)

% Generate matrix.
E_mRNA_degradation = MatrixGeneration;

% Generate a list of TUs to be degraded.
TUlist = Gene_list.TU;
TU_degraded = TUlist(~ismember(TUlist,TU_cleaved));
TU_degraded = cat(1,TU_degraded,Sub_TU_cleaved.ID);

% mRNA dilution, including TUs and Sub-TUs.
for i = 1:length(TU_degraded)
    name = TU_degraded{i};
    RxnID = strcat('R_dilution_mRNA_',name);
    Compo = {strcat(name,'_c')};
    Coeff = -1;
    Rev = 0;
    Subproc = 'RNA dilution';
    Catalyst = '';

    E_RNA_dilution = MatrixOneReactionAdding(E_RNA_dilution,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);    
end

% Add formation of degradosome with oligoribonuclease.
Orn = FindFunctionalProtein('llmg_1825',Functional_protein);
CshA = FindFunctionalProtein('llmg_0369',Functional_protein);
Eno = FindFunctionalProtein('llmg_0617',Functional_protein);
PfkA = FindFunctionalProtein('llmg_1118',Functional_protein);
Pnp = FindFunctionalProtein('llmg_2044',Functional_protein);
RNase_J1 = FindFunctionalProtein('llmg_0302',Functional_protein);
RNase_J2 = FindFunctionalProtein('llmg_0876',Functional_protein);
RNase_Y = FindFunctionalProtein('llmg_2156',Functional_protein);
Product = 'mRNA_Degradation_Complex_Enzyme_c';

RxnID = 'R_mRNA_Degradation_Complex_Enzyme';
Compo = {Orn;CshA;Eno;PfkA;Pnp;RNase_J1;RNase_J2;RNase_Y;Product};
Coeff = [-1;-1;-1;-1;-1;-1;-1;-1;1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                 RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

% Add protein degradation reaction.    
RxnID = 'R_mRNA_Degradation_Complex_Enzyme_degradation';
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
RxnID = 'R_dilution_mRNA_Degradation_Complex_Enzyme';
Compo = {Product};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

% Add mRNA degradation reactions.
for i = 1:length(TU_degraded)
    TUname = TU_degraded{i};
    if contains(TUname,'_sub_')%is Sub-TU
        TUseq = GetSequence(TUname,Sub_TU_cleaved.ID,Sub_TU_cleaved.sequence);
        %determine the first base
        X = TUseq(1);
        B = {'A';'G';'C';'T'};
        b = {'M_atp_c';'M_gtp_c';'M_ctp_c';'M_utp_c'};
        x = GetSequence(X,B,b);
        %determine coefficient
        bases = CountDNABase(TUseq(2:end));
        nbases = length(TUseq);
        natp = round(nbases/4);
        nadp = natp;
        npi = natp;
        nh = natp+nbases-1;
        nh2o = nh;

        %add degradation reactions
        RxnID = strcat('R_',TUname,'_degradation');
        Compo = {strcat(TUname,'_c');'M_atp_c';'M_h2o_c';...
                 'M_amp_c';'M_gmp_c';'M_cmp_c';'M_ump_c';...
                 'M_h_c';'M_pi_c';'M_adp_c';x};
        Coeff = [-1;-natp;-nh2o;...
                 transpose(bases);...
                 nh;npi;nadp;1];
        Rev = 0;
        Subproc = 'mRNA degradation';
        Catalyst = 'mRNA_Degradation_Complex_Enzyme_c';

        E_mRNA_degradation = MatrixOneReactionAdding(E_mRNA_degradation,...
                                         RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
    else%is not Sub-TU
        TUseq = GetSequence(TUname,TU.TUID,TU.sequence);
        %determine the first base
        X = TUseq(1);
        B = {'A';'G';'C';'T'};
        b = {'M_atp_c';'M_gtp_c';'M_ctp_c';'M_utp_c'};
        x = GetSequence(X,B,b);
        %determine coefficient
        bases = CountDNABase(TUseq(2:end));
        nbases = length(TUseq);
        natp = round(nbases/4);
        nadp = natp;
        npi = natp;
        nh = natp+nbases-1;
        nh2o = nh;
        
        %add degradation reactions        
        RxnID = strcat('R_',TUname,'_degradation');
        Compo = {strcat(TUname,'_c');'M_atp_c';'M_h2o_c';...
                 'M_amp_c';'M_gmp_c';'M_cmp_c';'M_ump_c';...
                 'M_h_c';'M_pi_c';'M_adp_c';x};
        Coeff = [-1;-natp;-nh2o;...
                 transpose(bases);...
                 nh;npi;nadp;1];
        Rev = 0;
        Subproc = 'mRNA degradation';
        Catalyst = 'mRNA_Degradation_Complex_Enzyme_c';

        E_mRNA_degradation = MatrixOneReactionAdding(E_mRNA_degradation,...
                                         RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
    end
end
