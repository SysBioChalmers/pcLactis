%% CheckInactiveEnzyme
%   This function checks if produced enzymes do not carry fluxes.
function [tf,enzyme_inact,f_enzyme_inact,f_inact_part] = CheckInactiveEnzyme(model,sol_full,factor_k)
% Input: model and flux distribution.
%        factor_k is the global saturation factor.
% Output: tf, a logical result. Logical 1 means at least one enzyme is
%             inactive or its part is inactive. Logical 0 means no inactive
%             enzymes.
%         enzyme_inact, list of such inactive enzymes.
%         f_enzyme_inact, total proportion of such inactive enzymes in
%                         biomass (g/g).
%         f_inact_part, inactive mass of inactive enzymes in biomass (g/g).

key_transporters = {'M_ALAt2r_1_fwd_Enzyme_c',...%ala in
                        'M_ALAt2r_2_fwd_Enzyme_c',...%ala in
                        'M_ARGt2r_Enzyme_c',...%arg in
                        'M_ARGORNt7_1_fwd_Enzyme_c',...%arg in
                        'M_ARGORNt7_2_fwd_Enzyme_c',...%arg in
                        'M_ARGabc_Enzyme_c',...%arg in
                        'M_ASNt2r_Enzyme_c',...%asn in
                        'M_ASPt2r_fwd_Enzyme_c',...%asp in
                        'M_CYSt2r_Enzyme_c',...%cys in
                        'M_GLUabc_1_Enzyme_c',...%glu in
                        'M_GLUabc_2_Enzyme_c',...%glu in
                        'M_GLUABUTt7_fwd_Enzyme_c',...%glu in
                        'M_GLNabc_1_Enzyme_c',...%gln in
                        'M_GLNabc_2_Enzyme_c',...%gln in
                        'M_GLYt2r_1_fwd_Enzyme_c',...%gly in
                        'M_GLYt2r_2_fwd_Enzyme_c',...%gly in
                        'M_HISt2r_Enzyme_c',...%his in
                        'M_ILEt2r_fwd_Enzyme_c',...%ile in
                        'M_LEUt2r_fwd_Enzyme_c',...%leu in
                        'M_LYSt2r_1_fwd_Enzyme_c',...%lys in
                        'M_LYSt2r_2_fwd_Enzyme_c',...%lys in
                        'M_METabc_1_Enzyme_c',...%met in
                        'M_METabc_2_Enzyme_c',...%met in
                        'M_METabc_3_Enzyme_c',...%met in
                        'M_METabc_4_Enzyme_c',...%met in
                        'M_METt2r_Enzyme_c',...%met in
                        'M_PHEt2r_fwd_Enzyme_c',...%phe in
                        'M_PROabc_Enzyme_c',...%pro in
                        'M_SERt2r_fwd_Enzyme_c',...%ser in
                        'M_THRt2r_fwd_Enzyme_c',...%thr in
                        'M_TRPt2r_fwd_Enzyme_c',...%trp in
                        'M_TYRt2r_fwd_Enzyme_c',...%tyr in
                        'M_VALt2r_fwd_Enzyme_c',...%val in
                        'M_PNTOabc_Enzyme_c',...%pantothenate in
                        'M_PNTOt2_fwd_Enzyme_c',...%pantothenate in
                        'M_NACt_Enzyme_c',...%nicotinate in
                        'M_THMabc_Enzyme_c',...%thiamin in
                        'M_H2Ot_fwd_Enzyme_c',...%h2o in
                        'M_H2Ot_rvs_Enzyme_c',...%h2o out
                        'M_PIabc_Enzyme_c',...%pi in
                        'M_NH4t_fwd_Enzyme_c',...%nh4 in
                        'M_SO4t2_Enzyme_c',...%so4 in
                        'M_FE2abc_Enzyme_c',...%fe2 in
                        'M_FE3abc_Enzyme_c',...%fe3 in
                        'M_H2St1_rvs_Enzyme_c',...%h2s in
                        'M_MNabc_Enzyme_c',...%mn2 in
                        'M_MNt2_1_Enzyme_c',...%mn2 in
                        'M_MNt2_2_Enzyme_c',...%mn2 in
                        'M_ZNabc_Enzyme_c',...%zn2 in
                        };

load('Info_enzyme.mat');

load('kcat_M.mat');
m_enzyme = kcat_M.gpr;
m_kcat = kcat_M.kcat;

mu = sol_full(strcmp(model.rxns,'R_biomass_dilution'));

% Dilution of enzymes.
idx = cellfun(@(x) isequal({'Enzyme dilution'},x),model.subSystems,...
              'UniformOutput',false);
enzyme_dil_rxns = model.rxns(cell2mat(idx));
enzyme_dil_fluxes = sol_full(cell2mat(idx));
enzyme_dil_nonzero_fluxes_rxns = enzyme_dil_rxns(enzyme_dil_fluxes ~= 0);

%remove glucose uptake reaction
enzyme_dil_nonzero_fluxes_rxns = enzyme_dil_nonzero_fluxes_rxns(~contains(enzyme_dil_nonzero_fluxes_rxns,'R_dilution_M_GLCpts_2_Enzyme'));
%only keep metabolic reactions
enzyme_dil_nonzero_fluxes_rxns = enzyme_dil_nonzero_fluxes_rxns(contains(enzyme_dil_nonzero_fluxes_rxns,'R_dilution_M_'));

tf_list = [];
enzyme_inact = {};
f_inact_part = [];

% There could be a cycle between forward and backward reaction to carry
% inactive enzyme
idx_tmp = contains(enzyme_dil_nonzero_fluxes_rxns,'_rvs_') | contains(enzyme_dil_nonzero_fluxes_rxns,'_fwd_');
reversiblerxnid = enzyme_dil_nonzero_fluxes_rxns(idx_tmp);
reversiblerxn = cellfun(@(x) x(14:end-11),reversiblerxnid,'UniformOutput',false);
[~, tmp, ~] = unique(reversiblerxn);
reversiblerxnlist = reversiblerxn(setdiff((1:length(reversiblerxn))',tmp));
for i = 1:length(reversiblerxnlist)
    id_tmp = reversiblerxnlist{i};
    dil_rxn_fwd = strcat('R_dilution_M_',id_tmp,'_fwd_Enzyme');
    enzyme_id_fwd = strcat('M_',id_tmp,'_fwd_Enzyme_c');
    dil_rxn_rvs = strcat('R_dilution_M_',id_tmp,'_rvs_Enzyme');
    enzyme_id_rvs = strcat('M_',id_tmp,'_rvs_Enzyme_c');
    mflux_fwd = sol_full(strcmp(model.grRules,enzyme_id_fwd));
    mflux_rvs = sol_full(strcmp(model.grRules,enzyme_id_rvs));
    dilflux_fwd = sol_full(strcmp(model.rxns,dil_rxn_fwd));
    dilflux_rvs = sol_full(strcmp(model.rxns,dil_rxn_rvs));
    if mflux_fwd ~= 0 && mflux_rvs ~= 0
        netflux = mflux_fwd - mflux_rvs;
        
        if netflux > 0
            rxn_id = dil_rxn_fwd;
            enzyme_id = enzyme_id_fwd;
            dil_flux = dilflux_fwd;
            enzyme_id_totally_inactive = enzyme_id_rvs;
            dil_flux_totally_inactive = dilflux_rvs;
        else
            rxn_id = dil_rxn_rvs;
            enzyme_id = enzyme_id_rvs;
            dil_flux = dilflux_rvs;
            enzyme_id_totally_inactive = enzyme_id_fwd;
            dil_flux_totally_inactive = dilflux_fwd;
        end
        
        kcat = m_kcat(ismember(m_enzyme,enzyme_id));
        % Change kcats extremely low or high value %%%%%%%%%%%%%%%%%%%%%%%%
        if kcat < quantile(m_kcat,0.1,1)
            kcat = quantile(m_kcat,0.1,1);
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ismember(enzyme_id,key_transporters)
            factor_k_tmp = 1;
        else
            factor_k_tmp = factor_k;
        end
        if factor_k_tmp < 0
            factor_k_tmp = 0.01;
        elseif factor_k_tmp > 1
            factor_k_tmp = 1;
        end
        kcat = kcat * factor_k_tmp;
        
        required_enzyme_conc = netflux/kcat;
        real_enzyme_conc = dil_flux/mu;
        delta_enzyme_conc = real_enzyme_conc - required_enzyme_conc;
        
        tf_i = delta_enzyme_conc > 0;
        tf_list = [tf_list;tf_i];
        enzyme_inact = [enzyme_inact;enzyme_id];
        MW = Info_enzyme.MW(contains(Info_enzyme.ID,enzyme_id));
        f_i = delta_enzyme_conc*MW/1000;
        f_inact_part = [f_inact_part;f_i];
        
        enzyme_conc_totally_inactive = dil_flux_totally_inactive/mu;
        tf_i = enzyme_conc_totally_inactive > 0;
        tf_list = [tf_list;tf_i];
        enzyme_inact = [enzyme_inact;enzyme_id_totally_inactive];
        MW = Info_enzyme.MW(contains(Info_enzyme.ID,enzyme_id_totally_inactive));
        f_i = enzyme_conc_totally_inactive*MW/1000;
        f_inact_part = [f_inact_part;f_i];
    end
end

for i = 1:length(enzyme_dil_nonzero_fluxes_rxns)
    rxn_id = enzyme_dil_nonzero_fluxes_rxns{i};
    enzyme_id = strcat(rxn_id(12:end),'_c');
    
    if ~ismember(enzyme_id,enzyme_inact)
    
        kcat = m_kcat(ismember(m_enzyme,enzyme_id));
        % Change kcats extremely low or high value %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if kcat < quantile(m_kcat,0.1,1)
            kcat = quantile(m_kcat,0.1,1);
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ismember(enzyme_id,key_transporters)
            factor_k_tmp = 1;
        else
            factor_k_tmp = factor_k;
        end
        if factor_k_tmp < 0
            factor_k_tmp = 0.01;
        elseif factor_k_tmp > 1
            factor_k_tmp = 1;
        end
        kcat = kcat * factor_k_tmp;

        enzyme_carry_fluxes = sol_full(strcmp(model.grRules,enzyme_id));
        used_enzyme_conc = enzyme_carry_fluxes/kcat;

        dil_flux = sol_full(strcmp(model.rxns,rxn_id));
        real_enzyme_conc = dil_flux/mu;

        delta_enzyme_conc = abs(real_enzyme_conc - used_enzyme_conc);

        tf_i = delta_enzyme_conc > 0;
        tf_list = [tf_list;tf_i];

        if tf_i
            enzyme_inact = [enzyme_inact;enzyme_id];
            MW = Info_enzyme.MW(contains(Info_enzyme.ID,enzyme_id));
            f_i = delta_enzyme_conc*MW/1000;
            f_inact_part = [f_inact_part;f_i];
        end
    end
end

tf = sum(tf_list) > 0;
f_enzyme_inact = sum(f_inact_part);
