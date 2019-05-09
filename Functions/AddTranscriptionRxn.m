%% AddTranscriptionRxn
%   This function will generate transcription reactions representing
%   transcription elongation and termination.
function [E_Transcription,...
          E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = AddTranscriptionRxn(E_Transcription,Gene_list,...
                                                       TU,Functional_protein,...
                                                       E_Enzyme_formation,...
                                                       E_Enzyme_dilution,...
                                                       E_Protein_degradation)

% Add formation of transcription elongation and termination factors.
NusA = FindFunctionalProtein('llmg_1796',Functional_protein);
NusB = FindFunctionalProtein('llmg_1878',Functional_protein);
NusG = FindFunctionalProtein('llmg_2388',Functional_protein);
GreA = FindFunctionalProtein('llmg_0610',Functional_protein);
Mfd = FindFunctionalProtein('llmg_0013',Functional_protein);
Product = 'Transcription_Complex_Enzyme_c';

RxnID = 'R_Transcription_Complex_Enzyme';
Compo = {NusA;NusB;NusG;GreA;Mfd;Product};
Coeff = [-1;-1;-1;-1;-1;1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                 RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

% Add protein degradation reaction.    
RxnID = 'R_Transcription_Complex_Enzyme_degradation';
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
RxnID = 'R_dilution_Transcription_Complex_Enzyme';
Compo = {Product};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

% Add transcription reactions.
TUlist = Gene_list.TU;

for i = 1:length(TUlist)
    TUname = TUlist{i};
    index_TU = strcmp(TUname,TU.TUID(:,1));
    seq = TU.sequence{index_TU};
    bases = CountDNABase(seq(17:end));
    nppi = length(seq(17:end));

    RxnID = strcat('R_',TUname);
    Compo = {'M_atp_c';'M_gtp_c';'M_ctp_c';'M_utp_c';
             strcat(TUname,'_Initiated_c');...
             'M_ppi_c';strcat(TUname,'_c')};
    Coeff = [-transpose(bases);-1;nppi;1];
    Rev = 0;
    Subproc = 'Transcription';
    Catalyst = 'Transcription_Complex_Enzyme_c';

    E_Transcription = MatrixOneReactionAdding(E_Transcription,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
end
