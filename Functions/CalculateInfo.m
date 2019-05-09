%% CalculateInfo
%   This function will calculate length and molecular weight of RNA,
%   protein and enzymes, and 70S ribosome composition. 
function [Info_mRNA,...
          Info_tRNA,...
          Info_protein,...
          Info_enzyme,...
          Info_ribosome] = CalculateInfo(TU,TU_degraded,Sub_TU_cleaved,...
                                         Gene_list,Seq_gene_NA_total,...
                                         Translated_gene,Functional_protein,...
                                         E_Enzyme_formation)

% mRNA.
Info_mRNA = struct();
Info_mRNA.ID = TU_degraded;
Info_mRNA.seq = cell(0,1);
Info_mRNA.len = zeros(0,1);
Info_mRNA.MW = zeros(0,1);

for i = 1:length(Info_mRNA.ID)
    tu = Info_mRNA.ID{i};
    if ~isempty(find(strcmp(TU.TUID,tu),1))
        index = strcmp(TU.TUID,tu);
        seq_DNA = TU.sequence{index};
	else
        index = strcmp(Sub_TU_cleaved.ID,tu);
        seq_DNA = Sub_TU_cleaved.sequence{index};
    end
	seq = strrep(seq_DNA,'T','U');
	Info_mRNA.seq{i,1} = seq;
	Info_mRNA.len(i,1) = length(seq);
	Info_mRNA.MW(i,1) = getRNAMW(seq);
end

% tRNA.
Info_tRNA = struct();
Info_tRNA.ID = Gene_list.E_tRNA_gene;
Info_tRNA.seq = cell(0,1);
Info_tRNA.len = zeros(0,1);
Info_tRNA.MW = zeros(0,1);

for i = 1:length(Info_tRNA.ID)
    trna = Info_tRNA.ID{i};
    index = strcmp(Seq_gene_NA_total.gene,trna);
    seq_DNA = Seq_gene_NA_total.sequence{index};
	seq = strrep(seq_DNA,'T','U');
	Info_tRNA.seq{i,1} = seq;
	Info_tRNA.len(i,1) = length(seq);
	Info_tRNA.MW(i,1) = getRNAMW(seq);
end

% Protein.
%   protein -> protein_mature -> protein_folded -> enzyme
Info_protein = struct();
Info_protein.ID = Translated_gene.ID;
Info_protein.seq = Translated_gene.AAseq;
Info_protein.seq_mature = Translated_gene.AAseq_mature;
Info_protein.len = zeros(0,1);
Info_protein.len_mature = zeros(0,1);
Info_protein.stoich = zeros(0,1);
Info_protein.MW_mature = zeros(0,1);
Info_protein.MW_folded = zeros(0,1);
Info_protein.ID_folded = cell(0,1);

for i = 1:length(Info_protein.ID)
    protein = Info_protein.ID{i};
    Info_protein.len(i,1) = length(Info_protein.seq{i,1});
    Info_protein.len_mature(i,1) = length(Info_protein.seq_mature{i,1});
    index = strcmp(Functional_protein(:,1),protein);
    prot_stoich = Functional_protein{index,2};
    Info_protein.ID_folded{i,1} = prot_stoich;
    if contains(prot_stoich,'Monomer')
        st = 1;
    else
        z = strfind(prot_stoich,'mer');
        a = max(strfind(prot_stoich(1:z),'_'));
        st = str2double(prot_stoich(a+1:z-1));
    end
    Info_protein.stoich(i,1) = st;
    seq_mature = Info_protein.seq_mature{i,1};
    Info_protein.MW_mature(i,1) = getProteinMW(seq_mature);
    Info_protein.MW_folded(i,1) = st * Info_protein.MW_mature(i,1);  
end

% Enzyme.
Info_enzyme = struct();
Info_enzyme.ID = cell(0,1);
Info_enzyme.subunit = cell(0,1);
Info_enzyme.MW = zeros(0,1);

rxns = unique(E_Enzyme_formation.RxnList);
for i = 1:length(rxns)
    rxnid = rxns{i};
    idx = strcmp(E_Enzyme_formation.RxnList,rxnid);
    reactant = E_Enzyme_formation.CompoList(idx);
    coeff = E_Enzyme_formation.CoeffList(idx);
    enzymeid = reactant(coeff>0);
    subunitid = reactant(coeff<0);
    MW = 0;
    for j = 1:length(subunitid)
        prot = subunitid{j};
        index_coeff = strcmp(reactant,prot);
        coeff_subunit = -1*coeff(index_coeff);
        coeff_subunit = sum(coeff_subunit);
        index = strcmp(Info_protein.ID_folded,prot);
        MW = MW + coeff_subunit*Info_protein.MW_folded(index);
    end
    Info_enzyme.ID{i,1} = cell2mat(enzymeid);
    Info_enzyme.subunit{i,1} = subunitid;
    Info_enzyme.MW(i,1) = MW;
end

% 70S ribosome.
Info_ribosome = struct();

cd Data;
[~, ~, RibLSUinfo] = xlsread('Ribosomal_protein.xlsx','LSU');
[~, ~, RibSSUinfo] = xlsread('Ribosomal_protein.xlsx','SSU');
cd ../;

Info_ribosome.subunit = [RibLSUinfo(2:end,1);RibSSUinfo(2:end,1)];
Info_ribosome.coeff = cell2mat([RibLSUinfo(2:end,3);RibSSUinfo(2:end,3)]);

MW_protein = 0;
for i = 1:length(Info_ribosome.subunit)
    prot = Info_ribosome.subunit{i};
    coef = Info_ribosome.coeff(i);
    index = strcmp(Info_protein.ID,prot);
    MW_protein = MW_protein + coef*Info_protein.MW_folded(index);
end
Info_ribosome.MW_protein = MW_protein;

%5S-rRNA
index_5S = strcmp(Seq_gene_NA_total.gene,'llmg_rRNA_48b');
seq_DNA_5S = Seq_gene_NA_total.sequence{index_5S};
seq_5S = strrep(seq_DNA_5S,'T','U');
MW_5S = getRNAMW(seq_5S);
%16S-rRNA
index_16S = strcmp(Seq_gene_NA_total.gene,'llmg_rRNA_48a');
seq_DNA_16S = Seq_gene_NA_total.sequence{index_16S};
seq_16S = strrep(seq_DNA_16S,'T','U');
MW_16S = getRNAMW(seq_16S);
%23S-rRNA
index_23S = strcmp(Seq_gene_NA_total.gene,'llmg_rRNA_48c');
seq_DNA_23S = Seq_gene_NA_total.sequence{index_23S};
seq_23S = strrep(seq_DNA_23S,'T','U');
MW_23S = getRNAMW(seq_23S);

MW_RNA = MW_5S + MW_16S + MW_23S;

Info_ribosome.MW_RNA = MW_RNA;

Info_ribosome.MW_protein = MW_protein + MW_RNA;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MW = getRNAMW(seq)
A = 329; U = 306; G = 345; C = 305;
nA = length(strfind(seq,'A'));
nU = length(strfind(seq,'U'));
nG = length(strfind(seq,'G'));
nC = length(strfind(seq,'C'));
MW = A*nA + U*nU + G*nG + C*nC + 159;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MW = getProteinMW(seq)
% (From GECKO)
aa_codes = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N', ...
            'O','P','Q','R','S','T','U','V','W','X','Y','Z'};
aa_MWs   = [71.08 114.60 103.14 115.09 129.11 147.17 57.05 137.14 ...
            113.16 113.16 128.17 113.16 131.20 114.10 255.31 97.12 ...
            128.13 156.19 87.08 101.10 150.04 99.13 186.21 126.50 ...
            163.17 128.62];
MW = 18;
for i = 1:length(aa_codes)
    count = length(strfind(seq,aa_codes{i}));
    MW = MW + count*aa_MWs(i);
end
end

