%% Collect membrane proteins

% As proteomics data contains membrane and soluble proteins, here we
% collect membrane proteins in the model, which should catalyze the
% exchange of metabolites between the membrane.

% The data will be saved in 'Membrane_proteins.mat'

load('pcLactis_Model.mat');
load('Info_enzyme.mat');
model = pcLactis_Model;
mem_proteins = {};
idx = cellfun(@(x) isequal({'Metabolism'},x),model.subSystems,'UniformOutput',false);
M_rxns_idx = find(cell2mat(idx));
for i = 1:length(M_rxns_idx)
    rxnidx = M_rxns_idx(i);
	metlist = model.mets(model.S(:,rxnidx) ~= 0);
	if find(cell2mat(cellfun(@(x) strcmp(x(length(x)-1:end),'_e'),metlist,'UniformOutput',false)),1) > 0
        enzyme_id = model.grRules{rxnidx};
        if ~isempty(enzyme_id)
            mem_proteins = [mem_proteins;Info_enzyme.subunit{contains(Info_enzyme.ID,enzyme_id)}];
        end
	end
end
mem_proteins = unique(mem_proteins);

for i = 1:length(mem_proteins)
    name = mem_proteins{i};
    idx = strfind(name,'_');
    idx = idx(2);
    name = name(1:idx-1);
    mem_proteins{i} = name;
end

save('mem_proteins.mat','mem_proteins');

clear;
