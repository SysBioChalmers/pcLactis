%% AssignMetabolicKcat
%   This function collect and assign kcat for metabolic enzymes.
function [RxnGPR, kcats] = AssignMetabolicKcat(RxnGPR,Metabolic_Matrix)

% EC numbers were downloaded from Uniprot, BioCyc and KEGG.

% Store EC number.
cd Data;
[~, ~, rawUniprot] = xlsread('EC_number.xlsx','Uniprot');
[~, ~, rawKEGG] = xlsread('EC_number.xlsx','KEGG');
cd ../;
EC = struct();
EC.Uniprot = rawUniprot(:,1:2);
EC.KEGG = rawKEGG(:,1:2);

% Import collected protein stoichiometry data.
cd Data;
fid = fopen('protein_stoichiometry_20181030.txt','r');
data = textscan(fid,'%s %f');
fclose(fid);
cd ../;
Protein_stoichiometry(:,1) = data{1};
Protein_stoichiometry(:,2) = num2cell(data{2});

% Use RxnGPR to retrieve RxnID, GeneID, EC number and Substrate.
%split complexes
RxnGPR.Gene = RxnGPR.GPR;
for i = 1:length(RxnGPR.GPR)
    z = RxnGPR.GPR{i};
    if ~isempty(strfind(z,'and'))%if it is a complex
        num = length(strfind(z,'and'))+1;
        RxnGPR.Gene(i,1:num) = strsplit(z,' and ');
    end
end
%search EC number
[a,b] = size(RxnGPR.Gene);
RxnGPR.EC = cell(a,b);
[x,y] = find(cellfun('isempty',RxnGPR.Gene) == 0);
for i = 1:length(x)
    z = RxnGPR.Gene{x(i),y(i)};
    if ismember(z,EC.Uniprot(:,1))%search EC in Uniprot
        index = find(strcmp(z,EC.Uniprot(:,1)));
        RxnGPR.EC{x(i),y(i)} = EC.Uniprot{index,2};
    else%if not found in Uniprot
        if ismember(z,EC.KEGG(:,1))%search EC in KEGG
            index = find(strcmp(z,EC.KEGG(:,1)));
            RxnGPR.EC{x(i),y(i)} = EC.KEGG{index,2};
        end
    end
end
%add substrate
%replace metID with full name
[~, ~, raw_metID] = xlsread('M_model.xlsx','metabolites');
MetName = raw_metID(2:end,1:2);
MetName(:,2) = cellfun(@(x) strrep(x,'_',' '),MetName(:,2),'UniformOutput',false);
n = length(Metabolic_Matrix.CompoList);
Metabolic_Matrix.MetName = cell(n,1);
for i = 1:n
    index = find(strcmp(Metabolic_Matrix.CompoList{i},MetName(:,1)));
    Metabolic_Matrix.MetName{i} = MetName{index,2};
end
%then add substrate
RxnGPR.Subs = cell(length(RxnGPR.Rxn),1);
for i = 1:length(RxnGPR.Rxn)
    z = find(strcmp(RxnGPR.Rxn{i},Metabolic_Matrix.RxnList(:,1)));
    z = z(find(Metabolic_Matrix.CoeffList(z) < 0));
    RxnGPR.Subs(i,1:length(z)) = transpose(Metabolic_Matrix.MetName(z));
end

% Retrieve kcat.
%kcat collection based on GECKO. Need Rxn_EC_Sub structure.
cd kcat_GECKO;

org_name = 'Lactococcus lactis subsp. cremoris MG1363';%type KEGG organism name
cd kcatDatabases;
load('PhylDist.mat');%load phylogenetic distance info
cd ../;
org_index = find_inKEGG(org_name,phylDistStruct.names);

[KCATcell, SAcell] = loadBRENDAdata;%BRENDA data for all the organisms
cd ../;

m_matrix = Metabolic_Matrix;%M_Matrix_2 should be generated before.

Rxns       = RxnGPR.Rxn;
substrates = RxnGPR.Subs;
EC_numbers = RxnGPR.EC;

[mM,nM]      = size(EC_numbers);
forw.kcats   = zeros(mM,nM);
forw.org_s   = zeros(mM,nM);
forw.rest_s  = zeros(mM,nM);
forw.org_ns  = zeros(mM,nM);
forw.rest_ns = zeros(mM,nM);
forw.org_sa  = zeros(mM,nM);
forw.rest_sa = zeros(mM,nM);
tot.queries  = 0;
tot.org_s    = 0;%match organism and substrate
tot.rest_s   = 0;%match the closest organism but match the substrate
tot.org_ns   = 0;%match organism but with any substrate
tot.rest_ns  = 0;%match any organism and any substrate
tot.org_sa   = 0;%match organism but for any substrate (SA*MW)
tot.rest_sa  = 0;%match any organism and SA*MW
tot.wc0      = 0;
tot.wc1      = 0;
tot.wc2      = 0;
tot.wc3      = 0;
tot.wc4      = 0;
tot.matrix   = zeros(6,5);
%Main loop: 
for i = 1:mM
%Match:
    for j = 1:nM
        EC = EC_numbers{i,j};
        if ~isempty(EC)
            [forw,tot] = iterativeMatch(EC,[],substrates(i,:),i,j,...
                                    KCATcell,forw,tot,m_matrix,org_name,...
                                    phylDistStruct,org_index,SAcell,Rxns);
        end
    end
     %Display progress:
     disp(['Matching kcats: Ready with rxn ' num2str(i)])
 end
kcats = forw;
kcats.tot = tot;

% kcat times stoichiometry
[a,b] = size(kcats.kcats);
kcat_xstoich = zeros(a,b);

for i = 1:a
    for j = 1:b
        if ~isempty(RxnGPR.Gene{i,j}) && ismember(RxnGPR.Gene(i,j),Protein_stoichiometry(:,1))
            index = strcmp(RxnGPR.Gene(i,j),Protein_stoichiometry(:,1));
            stoich = Protein_stoichiometry{index,2};
            kcat_xstoich(i,j) = kcats.kcats(i,j)*stoich;
        else
            kcat_xstoich(i,j) = kcats.kcats(i,j);
        end
    end
end


% Assign kcat.
% For complex, the median kcat among the subunits were chosen as the kcat.
kcat_process = zeros(0,1);
for i = 1:length(kcat_xstoich)
    k = kcat_xstoich(i,:);
    if sum(abs(k(:))) == 0
        kcat_process(i,1) = 0;
    else 
        kcat_process(i,1) = min(k(find(k)));
%         kcat_process(i,1) = median(k(find(k)));
    end
end
% For the reaction without assigned kcat, median value of the kcats of
% assigned reactions was used. 
median_kcat = median(kcat_process(find(kcat_process)));
RxnGPR.kcat = zeros(0,1);
for i = 1:length(kcat_process)
    k = kcat_process(i,1);
    if k == 0
        RxnGPR.kcat(i,1) = median_kcat;
    else
        RxnGPR.kcat(i,1) = k;
    end
end
cd ../;
%export to excel file
kcat_xlsx(:,1) = RxnGPR.Rxn;
kcat_xlsx(:,2) = RxnGPR.newGPR;
kcat_xlsx(:,3) = num2cell(RxnGPR.kcat);
kcat_xlsx = cat(1,{'MRxnID','MRxnEnzyme','kcat(1/h)'},kcat_xlsx);
xlwrite('kcat_M.xls',kcat_xlsx);
