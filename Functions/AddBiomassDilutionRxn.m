%% AddBiomassDilutionRxn
%   This function will add a biomass dilution reaction.
function Matrix = AddBiomassDilutionRxn(Matrix,f_unmodeled)

% Set parameters.
% f_unmodeled, unmodeled protein proportion (g/gCDW)
GAM = 42; %GAM (mmolATP/gCDW/h)

% Calculate the stoichiometric coefficient of unmodeled protein in the
% biomass dilution reaction.
MW_unmodeled = 28290.75; %(g/mol)
S_unmodeled = 0.46*f_unmodeled/MW_unmodeled*1000;

% Import biomass data.
cd Data;
[~, ~, biomass_raw] = xlsread('Biomass.xlsx');
cd ../;

% Add biomass dilution reaction.
RxnID = 'R_biomass_dilution';
Sub = cat(1,biomass_raw(2:end,1),{'unmodeled_protein_biomass_c';'M_h2o_c';'M_atp_c'});
Prod = {'M_pi_c';'M_adp_c';'M_h_c'};
Compo = cat(1,Sub,Prod);
Coeff = [-1*cell2mat(biomass_raw(2:end,2));-1*[S_unmodeled;GAM;GAM];
         GAM;GAM;GAM];
Rev = 0;
Subproc = '';
Catalyst = '';
Matrix = MatrixOneReactionAdding(Matrix,RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

