%% CBMpy2COBRA
function [M_Matrix_original, M_Model] = CBMpy2COBRA(num_rxn,txt_rxn,...
                                                    num_matrix,txt_matrix)

M_Model_info = struct();
M_Model_info.Rxn = txt_rxn(2:end,1);
M_Model_info.Rev = num_rxn(:,1);
M_Model_info.GPR = txt_rxn(2:end,6);

% Convert CBMpy model to COBRA model.
%	Change CBMpy structure to standard structure.
 %   CBMpy structure:
    %  R_23PDE2pp				
    %  substrates	1	M_23cump_e	1	M_h2o_e
    %  products	1	M_h_e	1	M_3ump_e
    %  R_23PDE4pp				
    %  substrates	1	M_23ccmp_e	1	M_h2o_e
    %  products	1	M_3cmp_e	1	M_h_e
 %   Standard structure:
    %  R_23PDE2pp	M_23cump_e	-1
    %  R_23PDE2pp	M_h2o_e	-1
    %  R_23PDE2pp	M_h_e	1
    %  R_23PDE2pp	M_3ump_e	1
    %  R_23PDE4pp	M_23ccmp_e	-1
    %  R_23PDE4pp	M_h2o_e	-1
    %  R_23PDE4pp	M_3cmp_e	1
    %  R_23PDE4pp	M_h_e	1

M_Matrix_original = MatrixGeneration;

for i = 1:length(M_Model_info.Rxn)
    RxnID = M_Model_info.Rxn{i};
    Rxn_name = strcat('R_M_',RxnID(3:end));
    [~, row_txt_RxnID] = ismember(RxnID,txt_matrix);
    row_txt_Sub = row_txt_RxnID+1;
    row_txt_Prod = row_txt_RxnID+2;
    row_num_Sub = row_txt_Sub-1;
    row_num_Prod = row_txt_Prod-1;
    Sub_list = txt_matrix(row_txt_Sub,3:end);
    Prod_list = txt_matrix(row_txt_Prod,3:end);
    Sub_coeff_list = num_matrix(row_num_Sub,:);
    Prod_coeff_list = num_matrix(row_num_Prod,:);
    
    Sub = transpose(Sub_list(~isnan(Sub_coeff_list)));
    Prod = transpose(Prod_list(~isnan(Prod_coeff_list)));
    Sub_coeff = -1*transpose(Sub_coeff_list(~isnan(Sub_coeff_list)));
    Prod_coeff = transpose(Prod_coeff_list(~isnan(Prod_coeff_list)));

    CompoID = cat(1,Sub,Prod);
    Coeff = [Sub_coeff;Prod_coeff];
    Rev = M_Model_info.Rev(i);
    Subproc = 'Metabolism';
    Catalyst = M_Model_info.GPR{i};
    
    M_Matrix_original = MatrixOneReactionAdding(M_Matrix_original,Rxn_name,CompoID,...
                                         Coeff,Rev,Subproc,Catalyst);
end

% Import to Excel file.
M_Model = Matrix2CBModel(M_Matrix_original);
% writeCbModel(M_Model,'xls','M_model_Original.xls');

% Add reversibility information to M_Model.
M_Model.rev = [];
for i = 1:length(M_Model.rxns)
    rxnid = M_Model.rxns{i};
    rxnid = strrep(rxnid,'_M_','_');
    [~, index] = ismember(rxnid,M_Model_info.Rxn);
    reversibility = M_Model_info.Rev(index);
    M_Model.rev = [M_Model.rev;reversibility];
end
