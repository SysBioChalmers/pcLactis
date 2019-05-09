%% Predict protein levels using predicted flux distributions

% Minimize glucose concentration

% The results have been saved in 'predicted_proteomics2.mat'.

% Timing: ~ 1300 s

tic;

%% Use the codes below to generate data. 
load('Cfd2_fluxes_without_sf.mat');
load('Cfd2_fluxes_with_sf.mat');

% load the model with the correct S matrix the solutions were generated.
load('pcLactis_Model.mat');
model = pcLactis_Model;
GAM = 42;
NGAM = 2.5;
f_unmodeled = 0.42;
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

save('predicted_proteomics2.mat','predicted');

toc;