%% CalculateProteomics 
function AbsProteins = CalculateProteomics(model,sol_full)
% The unit of concentration is mmol/gCDW, and the result shows the
% concentration of monomer.

idx = cellfun(@(x) isequal({'Enzyme dilution'},x),model.subSystems,...
              'UniformOutput',false);
enzy_dil_rxns = model.rxns(cell2mat(idx));

mu = sol_full(strcmp(model.rxns,'R_biomass_dilution'));

proteinlist = {};
conclist = [];

for i = 1:length(enzy_dil_rxns)
    dil_rxn_id = enzy_dil_rxns(i);
    dil_flux = sol_full(strcmp(model.rxns,dil_rxn_id));
    enzy_id = model.mets(model.S(:,strcmp(model.rxns,dil_rxn_id)) ~= 0);
    
    if ~strcmp(enzy_id,'ribosome_70S_c')
        form_rxn_id = model.rxns(model.S(strcmp(model.mets,enzy_id),:) > 0);
        sub_list = model.mets(model.S(:,strcmp(model.rxns,form_rxn_id)) < 0);
        for j = 1:length(sub_list)
            sub_enzy_id = sub_list{j};
            coeff = -model.S(strcmp(model.mets,sub_enzy_id),strcmp(model.rxns,form_rxn_id));
            if contains(sub_enzy_id,'Monomer')
                n_monomer = 1;
                sub_enzy_id_temp = strrep(sub_enzy_id,'_Monomer_c','');
                if contains(sub_enzy_id_temp,'assumed')
                    protein_name = strrep(sub_enzy_id_temp,'_assumed','');
                else
                    protein_name = sub_enzy_id_temp;
                end
            else
                temp_id = strrep(sub_enzy_id,'mer_c','');
                n_monomer = str2double(temp_id(max(strfind(temp_id,'_'))+1:end));
                protein_name = temp_id(1:max(strfind(temp_id,'_'))-1);
            end
            conc_sub_enzy = coeff*n_monomer*dil_flux/mu;
            disp(['Calculating protein concentration: ' protein_name]);
            proteinlist = [proteinlist;{protein_name}];
            conclist = [conclist;conc_sub_enzy];
        end
    else
        rxn_id_30S = {'R_ribosome_30S_protein'};
        sub_list_30S = model.mets(model.S(:,strcmp(model.rxns,rxn_id_30S)) < 0);
        for j = 1:length(sub_list_30S)
            sub_enzy_id = sub_list_30S{j};
            coeff = -model.S(strcmp(model.mets,sub_enzy_id),strcmp(model.rxns,form_rxn_id));
            if contains(sub_enzy_id,'Monomer')
                n_monomer = 1;
                sub_enzy_id_temp = strrep(sub_enzy_id,'_Monomer_c','');
                if contains(sub_enzy_id_temp,'assumed')
                    protein_name = strrep(sub_enzy_id_temp,'_assumed','');
                else
                    protein_name = sub_enzy_id_temp;
                end
            else
                temp_id = strrep(sub_enzy_id,'mer_c','');
                n_monomer = str2double(temp_id(max(strfind(temp_id,'_'))+1:end));
                protein_name = temp_id(1:max(strfind(temp_id,'_'))-1);
            end
            conc_sub_enzy = coeff*n_monomer*dil_flux/mu;
            disp(['Calculating protein concentration: ' protein_name]);
            proteinlist = [proteinlist;{protein_name}];
            conclist = [conclist;conc_sub_enzy];
        end
        rxn_id_50S = {'R_ribosome_50S_protein'};
        sub_list_50S = model.mets(model.S(:,strcmp(model.rxns,rxn_id_50S)) < 0);
        for j = 1:length(sub_list_50S)
            sub_enzy_id = sub_list_50S{j};
            coeff = -model.S(strcmp(model.mets,sub_enzy_id),strcmp(model.rxns,form_rxn_id));
            if contains(sub_enzy_id,'Monomer')
                n_monomer = 1;
                sub_enzy_id_temp = strrep(sub_enzy_id,'_Monomer_c','');
                if contains(sub_enzy_id_temp,'assumed')
                    protein_name = strrep(sub_enzy_id_temp,'_assumed','');
                else
                    protein_name = sub_enzy_id_temp;
                end
            else
                temp_id = strrep(sub_enzy_id,'mer_c','');
                n_monomer = str2double(temp_id(max(strfind(temp_id,'_'))+1:end));
                protein_name = temp_id(1:max(strfind(temp_id,'_'))-1);
            end
            conc_sub_enzy = coeff*n_monomer*dil_flux/mu;
            disp(['Calculating protein concentration: ' protein_name]);
            proteinlist = [proteinlist;{protein_name}];
            conclist = [conclist;conc_sub_enzy];
        end
    end
end

AbsProteins = struct();
AbsProteins.protein_id = unique(proteinlist);
AbsProteins.concentration = zeros(length(AbsProteins.protein_id),1);

for i = 1:length(AbsProteins.protein_id)
    proteinid = AbsProteins.protein_id(i);
    AbsProteins.concentration(i,1) = sum(conclist(strcmp(proteinid,proteinlist)));
end
