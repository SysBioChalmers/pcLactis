%% GetTranscriptionUnits
function [TU, Gene_list] = GetTranscriptionUnits(Gene_list)

[~, txtTUinfo, ~] = xlsread('Transcription_units.xlsx');


TU = struct();
TU.TUID = txtTUinfo(2:end,1);
TU.gene = txtTUinfo(2:end,2);
TU.sequence = txtTUinfo(2:end,3);

% Select TUs should be involved in the ME model.
%   Due to co-transcription, more genes will be transcribed compared to the
%   original genes, i.e., Gene_list.total_gene.
Gene_list.TU = cell(0,1);
for i = 1:length(Gene_list.total_gene)
    x = length(find(cellfun('isempty',Gene_list.TU) == 0));
    gene = Gene_list.total_gene{i};
    if ~strcmp(gene,'dummy')
        if ismember(gene,TU.gene)
            index = strcmp(TU.gene,gene);
            Gene_list.TU(x+1,1) = TU.TUID(index,1);
        else
            error(['The gene ',gene,' is not involved in any TUs.'])
        end
    end
end
Gene_list.TU = unique(Gene_list.TU);%count the number of TUs in ME model

% Count the number of total genes based on TUs.
Gene_list.total_gene_in_all_TUs = cell(0,1);
for i = 1:length(Gene_list.TU)
    x = length(find(cellfun('isempty',Gene_list.total_gene_in_all_TUs) == 0));
    TUname = Gene_list.TU{i};
    index = find(strcmp(TU.TUID,TUname));
    n = length(index);
    Gene_list.total_gene_in_all_TUs(x+1:x+n,1) = TU.gene(index,1);
end
Gene_list.total_gene_in_all_TUs = unique(Gene_list.total_gene_in_all_TUs);
