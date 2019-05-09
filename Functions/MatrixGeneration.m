%% MatrixGeneration
% Generate a blank matrix.
function Matrix = MatrixGeneration
Matrix = struct();
Matrix.RxnList = cell(0,1);
Matrix.CompoList = cell(0,1);
Matrix.CoeffList = zeros(0,1);
Matrix.RevList = zeros(0,1);
Matrix.SubprocList = cell(0,1);
Matrix.CatalystList = cell(0,1);