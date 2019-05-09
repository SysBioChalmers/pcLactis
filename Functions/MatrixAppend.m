%% MatrixAppend
%   This function can append matrix and will be used to integrate different
%   subprocesses.
function Matrix = MatrixAppend(Matrix,Matrix_append)

Matrix.RxnList = cat(1,Matrix.RxnList,Matrix_append.RxnList);
Matrix.CompoList = cat(1,Matrix.CompoList,Matrix_append.CompoList);
Matrix.CoeffList = [Matrix.CoeffList;Matrix_append.CoeffList];
Matrix.RevList = [Matrix.RevList;Matrix_append.RevList];
Matrix.SubprocList = cat(1,Matrix.SubprocList,Matrix_append.SubprocList);
Matrix.CatalystList = cat(1,Matrix.CatalystList,Matrix_append.CatalystList);
