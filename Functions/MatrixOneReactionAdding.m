%% MatrixOneReactionAdding
% Add information to a matrix.
function Matrix = MatrixOneReactionAdding(Matrix,RxnID,CompoID,Coeff,Rev,Subproc,Catalyst)
% Input: Matrix, a structure containing RxnList, CompoList, CoeffList,
%                RevList, SubprocList and CatalystList.
%        RxnID, a string of reaction ID.
%        CompoID, a list of component ID, in a column.
%        Coeff, coefficient of each component, in a column.
%        Rev, whether the reaction is reversible (1) or not (0).
%        Subproc, a string of subprocess name.
%        Catalyst, one functional protein catalyzing the reaction.
% Output: Matrix, an updated matrix.
x = length(find(cellfun('isempty',Matrix.RxnList) == 0));%count rows
for i = 1:length(CompoID)
    Matrix.RxnList(x+i,1) = {RxnID};
    Matrix.CompoList(x+i,1) = CompoID(i,1);
    Matrix.CoeffList(x+i,1) = Coeff(i,1);
    Matrix.RevList(x+i,1) = Rev;
    Matrix.SubprocList(x+i,1) = {Subproc};
    Matrix.CatalystList(x+i,1) = {Catalyst};
end