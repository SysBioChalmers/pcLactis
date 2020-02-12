%% Predict transcript levels using predicted flux distributions

% The results have been saved in 'predicted_transcriptomics.mat'.

% Timing: ~ 10 s

tic;

%% Use the codes below to generate data. 
load('Cfd_fluxes_without_sf.mat');
load('Cfd_fluxes_with_sf.mat');

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
predicted_transcriptomics_without_sf = [];
[~,n] = size(fluxes_global_saturation_factor_unchanged);
for i = 1:n
    sol_tmp = fluxes_global_saturation_factor_unchanged(:,i);
    % Calculate absolute transcript. (mmol/gCDW)
    AbsTranscripts = CalculateTranscriptomics(model,sol_tmp);
    predicted_transcriptomics_without_sf = [predicted_transcriptomics_without_sf,AbsTranscripts.concentration];
end
geneID = AbsTranscripts.gene_id;
clear n i sol_tmp AbsTranscripts;

% solutions with changed saturation factor.
predicted_transcriptomics_with_sf = [];
[~,n] = size(fluxes_global_saturation_factor_changed);
for i = 1:n
    sol_tmp = fluxes_global_saturation_factor_changed(:,i);
    % Calculate absolute transcript. (mmol/gCDW)
    AbsTranscripts = CalculateTranscriptomics(model,sol_tmp);
    predicted_transcriptomics_with_sf = [predicted_transcriptomics_with_sf,AbsTranscripts.concentration];
end
clear n i sol_tmp AbsTranscripts;

predicted = struct();
predicted.without_sf = predicted_transcriptomics_without_sf;
predicted.with_sf = predicted_transcriptomics_with_sf;
predicted.geneID = geneID;

save('predicted_transcriptomics.mat','predicted');

toc;