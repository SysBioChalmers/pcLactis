%% CheckComponentBalance
function [tf,unbalanced_compo,net_fluxes]= CheckComponentBalance(model,sol_full)
% Input: model and flux distribution.
% Output: tf, one means there is unbalanced component, zero means no;
%         unbalanced_compo.

compo_net_flux = zeros(length(model.mets),1);
for i = 1:length(model.mets)
    disp(['Check component balance: ' num2str(i) '/' num2str(length(model.mets))]);
    compo_id = model.mets{i};
	compo_net_flux(i) = CalculateMetBalance(model,model.rxns,sol_full,compo_id);
end

tf = find(compo_net_flux,1) > 0;
idx = compo_net_flux ~= 0;
unbalanced_compo = model.mets(idx);
net_fluxes = compo_net_flux(idx);
