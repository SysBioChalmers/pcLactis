%% CheckInactiveEnzyme
%   This function checks if produced enzymes do not carry fluxes.
function [tf,enzyme_inact,f_enzyme_inact] = CheckInactiveEnzyme(model,sol_full)
% Input: model and flux distribution.
% Output: tf, a logical result. Logical 1 means at least one enzyme is
%             produced but does not carry flux. Logical 0 means no such
%             enzymes.
%         enzyme_inact, list of such inactive enzymes.
%         f_enzyme_inact, total proportion of such inactive enzymes in
%                         biomass (g/g).

load('Info_enzyme.mat');

mu = sol_full(strcmp(model.rxns,'R_biomass_dilution'));

% Dilution of enzymes.
idx = cellfun(@(x) isequal({'Enzyme dilution'},x),model.subSystems,...
              'UniformOutput',false);
enzyme_dil_rxns = model.rxns(cell2mat(idx));
enzyme_dil_fluxes = sol_full(cell2mat(idx));
enzyme_dil_nonzero_fluxes_rxns = enzyme_dil_rxns(enzyme_dil_fluxes ~= 0);

tf_list = [];
enzyme_inact = {};
f_enzyme_inact = 0;

for i = 1:length(enzyme_dil_nonzero_fluxes_rxns)
    rxn_id = enzyme_dil_nonzero_fluxes_rxns{i};
    enzyme_id = strcat(rxn_id(12:end),'_c');
    enzyme_carry_fluxes = sol_full(strcmp(model.grRules,enzyme_id));
    
    tf_i = max(abs(enzyme_carry_fluxes)) == 0;
    tf_list = [tf_list;tf_i];
    
    if tf_i
        enzyme_inact = [enzyme_inact;enzyme_id];
        
        dil_flux = sol_full(strcmp(model.rxns,rxn_id));
        MW = Info_enzyme.MW(contains(Info_enzyme.ID,enzyme_id));
        f_i = dil_flux*MW/1000/mu;
        f_enzyme_inact = f_enzyme_inact + f_i;
    end
end

tf = sum(tf_list) > 0;
    
