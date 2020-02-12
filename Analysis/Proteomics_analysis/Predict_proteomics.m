%% Predict protein levels using predicted flux distributions

% The results have been saved in 'predicted_proteomics.mat'.

% Timing: ~ 1100 s

tic;

%% Use the codes below to generate data. 
load('Cfd_fluxes_with_sf.mat');
load('Cfd_fluxes_without_sf.mat');

% load the model with the correct S matrix the solutions were generated.
load('pcLactis_Model.mat');
model = pcLactis_Model;
GAM = 36; %ATP coefficient in the new biomass equation.
NGAM = 2; %(mmol/gCDW/h)
f_unmodeled = 0.4; %proportion of unmodeled protein in total protein (g/g)
[model,~] = ChangeUnmodeledProtein(model,f_unmodeled);
model = ChangeATPinBiomass(model,GAM);
model = changeRxnBounds(model,'R_M_ATPM',NGAM,'b');
clear f_unmodeled GAM pcLactis_Model NGAM;

% solutions with unchanged saturation factor.
predicted_proteomics_without_sf = [];
[~,n] = size(fluxes_global_saturation_factor_unchanged);
for i = 1:n
    sol_tmp = fluxes_global_saturation_factor_unchanged(:,i);
    % Calculate absolute protein. (mmol/gCDW)
    AbsProteins = CalculateProteomics(model,sol_tmp);
    predicted_proteomics_without_sf = [predicted_proteomics_without_sf,AbsProteins.concentration];
end
proteinID = AbsProteins.protein_id;
clear n i sol_tmp AbsProteins;

% solutions with changed saturation factor.
predicted_proteomics_with_sf = [];
[~,n] = size(fluxes_global_saturation_factor_changed);
for i = 1:n
    sol_tmp = fluxes_global_saturation_factor_changed(:,i);
    % Calculate absolute protein. (mmol/gCDW)
    AbsProteins = CalculateProteomics(model,sol_tmp);
    predicted_proteomics_with_sf = [predicted_proteomics_with_sf,AbsProteins.concentration];
end
clear n i sol_tmp AbsProteins;

predicted = struct();
predicted.without_sf = predicted_proteomics_without_sf;
predicted.with_sf = predicted_proteomics_with_sf;
predicted.proteinID = proteinID;

save('predicted_proteomics.mat','predicted');

toc;