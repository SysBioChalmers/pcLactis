%% AddTranscriptionInitiationRxn
function [E_Transcription,...
          E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation] = AddTranscriptionInitiationRxn(E_Transcription,...
                                                                 Gene_list,TU,...
                                                                 Seq_gene_NA_total,...
                                                                 Functional_protein,...
                                                                 E_Enzyme_formation,...
                                                                 E_Enzyme_dilution,...
                                                                 E_Protein_degradation)

% Check if all the TUs contain the intact sequence of each gene.
[result, ~] = CheckSequenceConsistency(Seq_gene_NA_total,TU);
if result == 0
    error('Gene sequence is not in line with TU sequence.');
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % By doing this, we found that for 'llmg_0478' NCBI and BioCyc have     %
 % different sequences, and then adopted the NCBI sequence in the model, %
 % and changed the corresponding TU sequence, i.e., 'TU1G1R-646'. And    %
 % also for 'llmg_2311'.                                                 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add formation of RNA polymerase with sigma factor.
alpha = FindFunctionalProtein('llmg_2354',Functional_protein);
beta = FindFunctionalProtein('llmg_1982',Functional_protein);
beta_prime = FindFunctionalProtein('llmg_1981',Functional_protein);
omega = FindFunctionalProtein('llmg_2154',Functional_protein);
delta = FindFunctionalProtein('llmg_0608',Functional_protein);
sigma_factor = FindFunctionalProtein('llmg_0521',Functional_protein);
RNAP_sf = 'RNAP_sf_Enzyme_c';

RxnID = 'R_RNAP_sf_Enzyme';
Compo = {alpha;beta;beta_prime;omega;delta;sigma_factor;RNAP_sf};
Coeff = [-2;-1;-1;-1;-1;-1;1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                 RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

% Add protein degradation reaction.    
RxnID = 'R_RNAP_sf_Enzyme_degradation';
Sub_raw = Compo(1:end-1);
Sub = cellfun(@(x) strcat(x(1:end-2),'_degradation_c'),...
                      Sub_raw,'UniformOutput',false);
Compo = cat(1,{RNAP_sf},Sub);
Coeff = [-1;-1*Coeff(1:end-1)];
Rev = 0;
Subproc = 'Protein degradation';
Catalyst = '';
E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

% Add dilution reaction.
RxnID = 'R_dilution_RNAP_sf_Enzyme';
Compo = {'RNAP_sf_Enzyme_c'};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                            
% Add transcription initiation reactions.
TUlist = Gene_list.TU;
for i = 1:length(TUlist)
    TUname = TUlist{i};
    index_TU = strcmp(TUname,TU.TUID(:,1));
    seq = TU.sequence{index_TU};
    bases = CountDNABase(seq(1:16));

    RxnID = strcat('R_',TUname,'_Initiation');
    Compo = {'M_atp_c';'M_gtp_c';'M_ctp_c';'M_utp_c';...
             'M_ppi_c';strcat(TUname,'_Initiated_c')};
    Coeff = [-transpose(bases);15;1];
    Rev = 0;
    Subproc = 'Transcription';
    Catalyst = 'RNAP_sf_Enzyme_c';

    E_Transcription = MatrixOneReactionAdding(E_Transcription,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
end
