%% AddRNACleaveRxn
%   This function will generate stable RNA cleave reactions for tRNA and
%   rRNA from TUs.
function [E_Stable_RNA_cleavage,...
          Sub_TU_cleaved,...
          E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation,...
          TU_cleaved] = AddRNACleaveRxn(Gene_list,TU,Seq_gene_NA_total,...
                                        Functional_protein,...
                                        E_Enzyme_formation,...
                                        E_Enzyme_dilution,...
                                        E_Protein_degradation)

% Generate matrix.
E_Stable_RNA_cleavage = MatrixGeneration;

% Add formation of RNase.
RNase = FindFunctionalProtein('llmg_1753',Functional_protein);
Product = 'RNase_Enzyme_c';
RxnID = 'R_RNase_Enzyme';
Compo = {RNase;Product};
Coeff = [-1;1];
Rev = 0;
Subproc = 'Enzyme formation';
Catalyst = '';

E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                 RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                             
% Add protein degradation reaction.    
RxnID = 'R_RNase_Enzyme_degradation';
Compo = {'RNase_Enzyme_c';strcat(RNase(1:end-2),'_degradation_c')};
Coeff = [-1;1];
Rev = 0;
Subproc = 'Protein degradation';
Catalyst = '';
E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

% Add dilution reaction.
RxnID = 'R_dilution_RNase_Enzyme';
Compo = {'RNase_Enzyme_c'};
Coeff = -1;
Rev = 0;
Subproc = 'Enzyme dilution';
Catalyst = '';
E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

% Add Stable RNA cleavage reactions for TUs containning tRNA or rRNA.
%	Determine stable RNA.
Stb_RNA = unique(cat(1,Gene_list.E_rRNA_gene,Gene_list.E_tRNA_gene));

%	Determine TUs to be cleaved.
TU_cleaved = {};
for i = 1:length(Stb_RNA)
    RNA = Stb_RNA(i,1);
    index = strcmp(RNA,TU.gene(:,1));
    TU_cleaved = cat(1,TU_cleaved,TU.TUID(index,1));
end
TU_cleaved = unique(TU_cleaved);

% Main part.
Sub_TU_cleaved = struct();
Sub_TU_cleaved.ID = {};
Sub_TU_cleaved.sequence = {};

n = length(TU_cleaved);
for i = 1:n
    TUname = TU_cleaved{i};
    TUseq = GetSequence(TUname,TU.TUID,TU.sequence);
    index_TU = find(strcmp(TUname,TU.TUID));
    m = length(index_TU);
    if m == 1%TU contains only one stable RNA gene
        genename = TU.gene{index_TU};
        %add cleavage reactions        
        RxnID = strcat('R_',TUname,'_cleavage');
        Compo = {strcat(TUname,'_c');strcat(genename,'_unmodified_c')};
        Coeff = [-1;1];
        Rev = 0;
        Subproc = 'Stable RNA cleavage';
        Catalyst = '';
        E_Stable_RNA_cleavage = MatrixOneReactionAdding(E_Stable_RNA_cleavage,...
                                         RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
    else%TU contains more than one genes
        genename = TU.gene(index_TU,1);
        index_stbgene = cell2mat((cellfun(@(x) ismember(x,Stb_RNA),...
                                       genename,'UniformOutput',false)));
        stbgene = genename(index_stbgene,1);%find stable RNA gene(s)
        index_gene = cell2mat((cellfun(@(x) ismember(x,stbgene),...
                                    Seq_gene_NA_total.gene,...
                                    'UniformOutput',false)));
        geneseq = Seq_gene_NA_total.sequence(index_gene,1);
        
        stbgene_inTU = {};
        subTU = {};    
        subTUseq = {};
        while ~isempty(TUseq)
            genepos = cell2mat((cellfun(@(x) strfind(TUseq,x),...
                                geneseq,'UniformOutput',false)));
            if min(genepos) > 1%sequence does not start with a stable RNA gene
                subTU = cat(1,subTU,{strcat(TUname,'_sub_',num2str(length(subTU)+1),'_c')});
                subTUseq = cat(1,subTUseq,TUseq(1:min(genepos)-1));
                TUseq = TUseq(min(genepos):end);
            elseif min(genepos) == 1%sequence starts with a stable RNA gene
                index_startgene = find(genepos == min(genepos));
                startgene = stbgene{index_startgene};
                startgeneseq = geneseq{index_startgene};
                stbgene = cat(1,stbgene(1:index_startgene-1,1),stbgene(index_startgene+1:end,1));
                geneseq = cat(1,geneseq(1:index_startgene-1,1),geneseq(index_startgene+1:end,1));
                stbgene_inTU = cat(1,stbgene_inTU,{strcat(startgene,'_unmodified_c')});
                nextpos = length(startgeneseq)+1;
                TUseq = TUseq(nextpos:end);
            elseif isempty(genepos)%sequence does not contain stable RNA genes
                subTU = cat(1,subTU,{strcat(TUname,'_sub_',num2str(length(subTU)+1),'_c')});
                subTUseq = cat(1,subTUseq,{TUseq});
                TUseq = '';
            end
        end
        
        Sub_TU_cleaved.ID = cat(1,Sub_TU_cleaved.ID,...
                                cellfun(@(x) x(1:end-2),subTU,...
                                        'UniformOutput',false));
        Sub_TU_cleaved.sequence = cat(1,Sub_TU_cleaved.sequence,subTUseq);
        
        %add cleavage reactions        
        RxnID = strcat('R_',TUname,'_cleavage');
        Compo = cat(1,{strcat(TUname,'_c')},{'M_h2o_c'},{'M_h_c'},stbgene_inTU,subTU);
        r = length(stbgene_inTU)+length(subTU);
        Coeff = [-1;-1*(r-1),;r-1;ones(r,1)];%add H2O and H
        Rev = 0;
        Subproc = 'Stable RNA cleavage';
        Catalyst = 'RNase_Enzyme_c';
        E_Stable_RNA_cleavage = MatrixOneReactionAdding(E_Stable_RNA_cleavage,...
                                         RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
    end
end

% Import to Excel file.
%model = Matrix2CBModel(E_Stable_RNA_cleavage);
%writeCbModel(model,'xls','Stable_RNA_cleavage.xls');     
