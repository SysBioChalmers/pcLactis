%% ChangeATPinBiomassCluster
%   This function will change GAM in the biomass dilution reaction.
function model = ChangeATPinBiomassCluster(model,GAM)

rxn_id = strcmp(model.rxns,'R_biomass_dilution');
% rxn_id = strcmp(model.rxns,'R_M_biomass_LLA');

h2o_id = strcmp(model.mets,'M_h2o_c');
atp_id = strcmp(model.mets,'M_atp_c');
pi_id = strcmp(model.mets,'M_pi_c');
adp_id = strcmp(model.mets,'M_adp_c');
h_id = strcmp(model.mets,'M_h_c');

model.S(h2o_id,rxn_id) = -GAM;
model.S(atp_id,rxn_id) = -GAM;
model.S(pi_id,rxn_id) = GAM;
model.S(adp_id,rxn_id) = GAM;
model.S(h_id,rxn_id) = GAM;