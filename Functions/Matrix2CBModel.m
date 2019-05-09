%% Matrix2CBModel
%   Convert a matrix to a cobra model.
function model = Matrix2CBModel(Matrix)
% Input: Matrix, a structure containing RxnList, CompoList, CoeffList,
%                RevList, SubprocList and CatalystList.
%        RxnID, a string of reaction ID.
%        CompoID, a list of component ID, in a column.
%        Coeff, coefficient of each component, in a column.
%        Rev, whether the reaction is reversible (1) or not (0).
%        Subproc, a string of subprocess name.
%        Catalyst, one functional protein catalyzing the reaction.
% Output: model, a cobra model.
model = struct();
model.rxns = cell(0,1);
model.lb = zeros(0,1);
model.ub = zeros(0,1);
model.mets = cell(0,1);
model.S = sparse(0,0);
model.b = zeros(0,1);
model.osense = -1;
model.rxnGeneMat = sparse(0,0);
model.c = zeros(0,1);
model.rules = cell(0,1);
model.genes = cell(0,1);
model.csense = char();

RxnID = unique(Matrix.RxnList);
%Using addReaction in COBRA to generate the model.
for i = 1:length(RxnID)
    j = find(strcmp(Matrix.RxnList,RxnID{i}));
    k = min(j);
    Compo = transpose(Matrix.CompoList(j));
    Coeff = transpose(Matrix.CoeffList(j));
    Rev = Matrix.RevList(k);
    Subproc = Matrix.SubprocList{k};
    Catalyst = Matrix.CatalystList{k};
    model = addReaction(model,...
    RxnID{i},...
    'metaboliteList',Compo,...
    'stoichCoeffList',Coeff,...
    'reversible',Rev,...
    'geneRule',Catalyst,...
    'subSystem',Subproc);
end
