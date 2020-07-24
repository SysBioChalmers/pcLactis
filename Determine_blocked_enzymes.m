% Determine blocked enzymes in M model

load('M_Model.mat');
model = M_Model;

blockedRxns = findBlockedReaction(model);
blockedRxns = blockedRxns';
blockedgrRules = model.grRules(ismember(model.rxns,blockedRxns));

otherRxns = model.rxns(~ismember(model.rxns,blockedRxns));
othergrRules = model.grRules(~ismember(model.rxns,blockedRxns));
otherEnzymes = cell(0,1);
for i = 1:length(othergrRules)
    gr = othergrRules(i);
    gr = cell2mat(gr);
    gr = strrep(gr,'(',''); gr = strrep(gr,')',''); 
    gr = strrep(gr,' or ',' '); gr = strrep(gr,' and ',' ');
    gr = strtrim(gr); gr = strsplit(gr);
    gr = gr';
    otherEnzymes = cat(1,otherEnzymes,gr);
end
otherEnzymes = unique(otherEnzymes);

blocked_enzymes = cell(0,1);
for i = 1:length(blockedgrRules)
    gr = blockedgrRules(i);
    gr = cell2mat(gr);
    gr = strrep(gr,'(',''); gr = strrep(gr,')',''); 
    gr = strrep(gr,' or ',' '); gr = strrep(gr,' and ',' ');
    gr = strtrim(gr); gr = strsplit(gr);
    
    for j = 1:length(gr)
        gene_tmp = gr(j);
        if ~ismember(gene_tmp,otherEnzymes)
            blocked_enzymes = cat(1,blocked_enzymes,gene_tmp);
        end
    end
end
blocked_enzymes = unique(blocked_enzymes);

Unblocked_dummy_proteins = sum(ismember(othergrRules,{'dummy'}))
