%% WriteLPSatFactorTmp

% This is just for reduced cost analysis of AA uptake.

% This is a new version of the WriteLP function. In this new function, the
% specific saturation factor of glucose transporter will be changed alone, 
% which is based on glucose concentration. The saturation factors of most 
% of the other enzymes are determined based on proteomics data if they 
% change linearly with growth rate. For the enzymes showing poor 
% correlation of saturation factor and growth rate (R < 0.8), stored in 
% sat_factor.mat, their saturation factors are assumed to be always 1. At 
% last, the saturation factors of the transporters of the unlimited 
% compounds in the medium are always set to be 1.

% Note that the unknown enzymes can be assumed to be unchanged or changed
% with growth rate. Line 151 152 182 183.


function fileName = WriteLPSatFactorTmp(model,mu,f,osenseStr,rxnID,factor_k,...
                                        f_transporter,kcat_glc,factor_glc,...
                                        Info_enzyme,...
                                        Info_mRNA,...
                                        Info_protein,...
                                        Info_ribosome,...
                                        Info_tRNA,...
                                        mu4parameter)
% f_transporter is proportion of glucose transporter in total proteome (Upper bound)

fileName = sprintf('Simulation.lp');
fptr = fopen(fileName,'w');

% Set objective function
osenseStr = strcat(osenseStr,'\n');
fprintf(fptr,osenseStr);
index_obj = find(ismember(model.rxns,rxnID));
fprintf(fptr,'obj: X%d\n',index_obj);
fprintf(fptr,'Subject To\n');

% Kinetic parameters
k = SetParameters(mu4parameter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) SV = 0.
% (From Ibrahim)
for i = 1:numel(model.mets)
    j = find(full(model.S(i,:)));
    for m = 1:numel(j)
        s = full(model.S(i,j(m)));
        if mod(m,200) == 0
            sep = newline;
        else
            sep = '';
        end
        if m == 1
           eq = sprintf('%.15f X%d',s,j(m));
        else
           if s>0
               eq = sprintf('%s + %.15f X%d%c',eq,s,j(m),sep);
           else
               eq = sprintf('%s %.15f X%d%c',eq,s,j(m),sep);
           end
        end
    end
    fprintf(fptr,'C%d: %s = 0\n',i,eq);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Coupling metabolic reactions and enzymes.

[~, ~, k_raw_m] = xlsread('k_parameter.xlsx','M');
[~, ~, k_raw_e] = xlsread('k_parameter.xlsx','E');

m_enzyme = k_raw_m(2:end,1);
m_kcat = cell2mat(k_raw_m(2:end,2));
e_enzyme = k_raw_e(2:end,1);
e_kcat = cell2mat(k_raw_e(2:end,2));

load('sat_factor.mat');

for i = 1:length(m_enzyme)
    enzyme = m_enzyme{i};
  	kcat = m_kcat(i);
    
% Set kcat for glucose transporter
    glc_transporter = {'M_GLCpts_1_Enzyme_c','M_GLCpts_2_Enzyme_c',...
                       'M_GLCt2_fwd_Enzyme_c','M_GLCt2_rvs_Enzyme_c'};
    if ismember(enzyme,glc_transporter)
        kcat = kcat_glc*3600;
    end
    
% Change kcats extremely low or high value
	if kcat < quantile(m_kcat,0.1,1)
        kcat = quantile(m_kcat,0.1,1);
	end
    
% Filter some transporters
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

	if ismember(enzyme,key_transporters) || ismember(enzyme,sat_factor.poor_Enzyme)
        factor_k_tmp = 1;
	elseif strcmp(enzyme,'M_GLCpts_1_Enzyme_c') || strcmp(enzyme,'M_GLCpts_2_Enzyme_c') || strcmp(enzyme,'M_GLCt2_fwd_Enzyme_c')
        factor_k_tmp = factor_glc;
    elseif ismember(enzyme,sat_factor.EnzymeID_pc)
        factor_k_tmp = factor_k;
    else
        factor_k_tmp = factor_k; %%%%%%%% The others to be changed
%         factor_k_tmp = 1; %%%%%%%% The others to be 1
	end
    if factor_k_tmp < 0
        factor_k_tmp = 0.01;
    elseif factor_k_tmp > 1
        factor_k_tmp = 1;
    end
    kcat = kcat * factor_k_tmp;
    
    %find enzyme formation reaction id
    id_syn = strcat('R_',enzyme(1:end-2));
    idx_syn = find(strcmp(model.rxns,id_syn));
    
    %calculate the coefficient of enzyme formation reaction rate
    coef = kcat/(k.deg_enzyme + mu);
    
    %find enzyme used reaction id (metabolic reaction)
    idx_rxn = find(ismember(model.grRules,enzyme));

    fprintf(fptr,'CM%d: X%d - %.15f X%d <= 0\n',...
                 i,idx_rxn,coef,idx_syn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Coupling some gene expression reactions and enzymes.

for i = 1:length(e_enzyme)
    enzyme = e_enzyme{i};
    %kcat = e_kcat(i);
%     kcat = 360000*factor_k; %%%%%%%% The others to be changed
    kcat = 360000; %%%%%%%% The others to be 1
    id_syn = strcat('R_',enzyme(1:end-2));
    idx_syn = find(strcmp(model.rxns,id_syn));
    
    coef = kcat/(k.deg_enzyme + mu);    

    idx_rxn = find(ismember(model.grRules,enzyme));
    
    %the enzyme catalyses one or more reactions
    if length(idx_rxn) == 1
%         fprintf(fptr,'CE%d: X%d - %.15f X%d <= 0\n',...
%                      i,idx_rxn,coef,idx_syn);
        fprintf(fptr,'CE%d: X%d - %.15f X%d = 0\n',...
                     i,idx_rxn,coef,idx_syn);
    elseif length(idx_rxn) > 1
        for j = 1
            eq = sprintf('X%d',idx_rxn(j));
        end
        for j = 2:length(idx_rxn)
            if mod(j,200) == 0
                sep = newline;
            else
                sep = '';
            end
            eq = sprintf('%s + X%d%c',eq,idx_rxn(j),sep);
        end
%         fprintf(fptr,'CE%d: %s - %.15f X%d <= 0\n',...
%                      i,eq,coef,idx_syn);
        fprintf(fptr,'CE%d: %s - %.15f X%d = 0\n',...
                     i,eq,coef,idx_syn);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) Coupling transcription reactions and RNA polymerase.

enzyme = 'RNAP_sf_Enzyme_c';
id_syn = 'R_RNAP_sf_Enzyme';
idx_syn = find(strcmp(model.rxns,id_syn));

kcat = k.cplx_rnap;
coef = kcat/(k.deg_enzyme + mu);

idx_rxn = find(ismember(model.grRules,enzyme));

for i = 1:length(idx_rxn)
    idx = idx_rxn(i);
    rxn_name = model.rxns{idx};
    tu_name = rxn_name(3:end-11);
    idx_tu = strcmp(Info_mRNA.ID,tu_name);
    len_tu = Info_mRNA.len(idx_tu);

	if mod(i,200) == 0
        sep = newline;
	else
        sep = '';
	end
    
    if i == 1
        eq = sprintf('%.15f X%d',len_tu,idx);
    else
        eq = sprintf('%s + %.15f X%d%c',eq,len_tu,idx,sep);
    end
end
% fprintf(fptr,'CRNAP: %s - %.15f X%d <= 0\n',eq,coef,idx_syn);
fprintf(fptr,'CRNAP: %s - %.15f X%d = 0\n',eq,coef,idx_syn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5) Coupling mRNA degradation reactions and mRNA degradation complex.

enzyme = 'mRNA_Degradation_Complex_Enzyme_c';
id_syn = 'R_mRNA_Degradation_Complex_Enzyme';
idx_syn = find(strcmp(model.rxns,id_syn));

kcat = k.cplx_mrnadeg;
coef = kcat/(k.deg_enzyme + mu);

idx_rxn = find(ismember(model.grRules,enzyme));

for i = 1:length(idx_rxn)
    idx = idx_rxn(i);
    rxn_name = model.rxns{idx};
    tu_name = rxn_name(3:end-12);
    idx_tu = strcmp(Info_mRNA.ID,tu_name);
    len_tu = Info_mRNA.len(idx_tu);

	if mod(i,200) == 0
        sep = newline;
	else
        sep = '';
	end
    
    if i == 1
        eq = sprintf('%.15f X%d',len_tu,idx);
    else
        eq = sprintf('%s + %.15f X%d%c',eq,len_tu,idx,sep);
    end
end
% fprintf(fptr,'CmRNAcplx: %s - %.15f X%d <= 0\n',eq,coef,idx_syn);
fprintf(fptr,'CmRNAcplx: %s - %.15f X%d = 0\n',eq,coef,idx_syn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6) Coupling tRNA.
% tRNA charging reactions should be coupled with both synthetase formation
% reactions and generic tRNA renaming reactions. The former has been
% added in 3), so we just add the latter below.

% Collect generic tRNA.
idx = cellfun(@(x) isequal({'Generic RNA renaming'},x),...
              model.subSystems,'UniformOutput',false);
Renaming_rxns = model.rxns(cell2mat(idx));
tRNA_rxns = Renaming_rxns(contains(Renaming_rxns,'to_tRNA'));
ge_tRNA = cellfun(@(x) x(strfind(x,'to_')+3:end),...
              tRNA_rxns,'UniformOutput',false);
ge_tRNA = unique(ge_tRNA);

kcat = k.trna;
coef = kcat/mu;

for i = 1:length(ge_tRNA)
    trna = ge_tRNA{i};
    idx_renaming = find(contains(model.rxns,strcat('_to_',trna)));
    idx_charging = find(contains(model.rxns,strcat('_charging_',trna)));
    
    if length(idx_renaming) == 1
        fprintf(fptr,'CtRNA%d: X%d - %.15f X%d <= 0\n',...
                     i,idx_charging,coef,idx_renaming);
%         fprintf(fptr,'CtRNA%d: X%d - %.15f X%d = 0\n',...
%                      i,idx_charging,coef,idx_renaming);
    else
        for j = 1:length(idx_renaming)
            idx = idx_renaming(j);
        if j == 1
            eq = sprintf('X%d - %.15f X%d',idx_charging,coef,idx);
        else
            eq = sprintf('%s - %.15f X%d',eq,coef,idx);
        end
        end
        fprintf(fptr,'CtRNA%d: %s <= 0\n',i,eq);
%         fprintf(fptr,'CtRNA%d: %s = 0\n',i,eq);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7) Coupling mRNA.
% Collect TUs for translation.
tu = Info_mRNA.ID(~contains(Info_mRNA.ID,'sub'));
idx = cellfun(@(x) isequal({'Translation'},x),model.subSystems,...
              'UniformOutput',false);
transl_id = model.rxns(cell2mat(idx));
tu_transl = {};
for i = 1:length(tu)
    tu_id = strcat(tu{i},'_');
    if ~isempty(find(contains(transl_id,tu_id), 1))
        tu_transl = [tu_transl;{tu_id}];
    end
end
tu_transl = cellfun(@(x) x(1:end-1),tu_transl,'UniformOutput',false);

kcat = k.mrna;
coef = kcat/(k.deg_mrna + mu);

for i = 1:length(tu_transl)
    tu = tu_transl{i};
    
    transl_init_id = strcat('R_translation_initiation_',tu,'_');
    idx_transl_init = find(contains(model.rxns,transl_init_id));
    
    transcrip_id = strcat('R_',tu);
    idx_transcrip = find(strcmp(model.rxns,transcrip_id));
    
	if length(idx_transl_init) == 1
        fprintf(fptr,'CmRNA%d: X%d - %.15f X%d <= 0\n',...
                     i,idx_transl_init,coef,idx_transcrip);
%         fprintf(fptr,'CmRNA%d: X%d - %.15f X%d = 0\n',...
%                      i,idx_transl_init,coef,idx_transcrip);
    else
        for j = 1:length(idx_transl_init)
            idx = idx_transl_init(j);
        if j == 1
            eq = sprintf('X%d',idx);
        else
            eq = sprintf('%s + X%d',eq,idx);
        end
        end
        fprintf(fptr,'CmRNA%d: %s - %.15f X%d <= 0\n',...
                     i,eq,coef,idx_transcrip);
%         fprintf(fptr,'CmRNA%d: %s - %.15f X%d = 0\n',...
%                      i,eq,coef,idx_transcrip);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8) Coupling translation reactions and ribosome.

enzyme = 'ribosome_70S_c';
id_syn = 'R_ribosome_70S';
idx_syn = find(strcmp(model.rxns,id_syn));

kcat = k.ribo;
coef = kcat/(k.deg_ribosome + mu);

id_rxn = model.rxns(ismember(model.grRules,enzyme));

for i = 1:length(Info_protein.ID)
    prot_name = Info_protein.ID{i};
    rxn_name = id_rxn(contains(id_rxn,prot_name));
    idx_rxn = find(ismember(model.rxns,rxn_name));
    prot_len = Info_protein.len(i);

	if mod(i,200) == 0
        sep = newline;
	else
        sep = '';
	end
    
    if i == 1
        eq = sprintf('%.15f X%d',prot_len,idx_rxn);
    else
        eq = sprintf('%s + %.15f X%d%c',eq,prot_len,idx_rxn,sep);
    end
end
fprintf(fptr,'Cribosome: %s - %.15f X%d <= 0\n',eq,coef,idx_syn);
% fprintf(fptr,'Cribosome: %s - %.15f X%d = 0\n',eq,coef,idx_syn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9) Coupling protein degradation reactions and protease.

enzyme = 'Clp_Protease_Complex_Enzyme_c';
id_syn = 'R_Clp_Protease_Complex_Enzyme';
idx_syn = find(strcmp(model.rxns,id_syn));

kcat = k.protease;
coef = kcat/(k.deg_enzyme + mu);

idx_rxn = find(ismember(model.grRules,enzyme));

for i = 1:length(idx_rxn)
    idx = idx_rxn(i);
    rxn_name = model.rxns{idx};
    prot_name = strrep(rxn_name,'R_','');
    prot_name = strrep(prot_name,'_degradation','');
    prot_name = prot_name(1:max(strfind(prot_name,'_'))-1);

    idx_prot = strcmp(Info_protein.ID,prot_name);
    len_prot = Info_protein.len_mature(idx_prot);

	if mod(i,200) == 0
        sep = newline;
	else
        sep = '';
	end
    
    if i == 1
        eq = sprintf('%.15f X%d',len_prot,idx);
    else
        eq = sprintf('%s + %.15f X%d%c',eq,len_prot,idx,sep);
    end
end
% fprintf(fptr,'Cprotease: %s - %.15f X%d <= 0\n',eq,coef,idx_syn);
fprintf(fptr,'Cprotease: %s - %.15f X%d = 0\n',eq,coef,idx_syn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10) Degradation of mRNA, 70S ribosome and enzymes.

% Degradation of mRNA.
idx = cellfun(@(x) isequal({'mRNA degradation'},x),model.subSystems,...
              'UniformOutput',false);
rna_deg_rxns = model.rxns(cell2mat(idx));

coef = k.deg_mrna/(k.deg_mrna + mu);

for i = 1:length(rna_deg_rxns)
    rxn_id = rna_deg_rxns{i};
    idx_deg = find(strcmp(model.rxns,rxn_id));
    if contains(rxn_id,'sub')
        rna_name = rxn_id(3:end);
        rna_name = strrep(rna_name,'_degradation','');
        rna_met_name = strcat(rna_name,'_c');
        met_idx = strcmp(model.mets,rna_met_name);
        idx_cleavage = find(full(model.S(met_idx,:)) > 0);        

        if length(idx_cleavage) == 1
            fprintf(fptr,'CRNAdeg%d: X%d - %.15f X%d = 0\n',...
                         i,idx_deg,coef,idx_cleavage);
        else
            for j = 1:length(idx_cleavage)
                idx = idx_cleavage(j);
            if j == 1
                eq = sprintf('X%d - %.15f X%d',idx_deg,coef,idx);
            else
                eq = sprintf('%s - %.15f X%d',eq,coef,idx);
            end
            end
            fprintf(fptr,'CRNAdeg%d: %s = 0\n',i,eq);
        end
    else
        rna_name = rxn_id(3:end);
        rna_name = strrep(rna_name,'_degradation','');
        transcription_id = strcat('R_',rna_name);
        idx_transcription = find(strcmp(model.rxns,transcription_id));        
        fprintf(fptr,'CRNAdeg%d: X%d - %.15f X%d = 0\n',...
                     i,idx_deg,coef,idx_transcription);
    end
end

% Degradation of enzymes.
enzyme_deg_rxns = model.rxns(contains(model.rxns,'_Enzyme_degradation'));

coef = k.deg_enzyme/(k.deg_enzyme + mu);

for i = 1:length(enzyme_deg_rxns)
    rxn_id = enzyme_deg_rxns{i};
    idx_deg = find(strcmp(model.rxns,rxn_id));
    enzyme_formation_id = strrep(rxn_id,'_degradation','');
    idx_formation = find(strcmp(model.rxns,enzyme_formation_id));
    fprintf(fptr,'Cenzymedeg%d: X%d - %.15f X%d = 0\n',...
                 i,idx_deg,coef,idx_formation);
end
% Degradation of unmodeled protein.
ump_deg_idx = find(strcmp(model.rxns,'R_unmodeled_protein_degradation'));
ump_syn_idx = find(strcmp(model.rxns,'R_unmodeled_protein'));
coef = k.deg_enzyme/(k.deg_enzyme + mu);
fprintf(fptr,'Cumpdeg: X%d - %.15f X%d = 0\n',...
              ump_deg_idx,coef,ump_syn_idx);

% Degradation of 70S ribosome.

coef = k.deg_ribosome/(k.deg_ribosome + mu);

ribosome_deg_rxn = 'R_ribosome_70S_degradation';
idx_deg = find(strcmp(model.rxns,ribosome_deg_rxn));
ribosome_formation_rxn = 'R_ribosome_70S';
idx_ribosome_formation = find(strcmp(model.rxns,ribosome_formation_rxn));
fprintf(fptr,'Cribodeg: X%d - %.15f X%d = 0\n',...
             idx_deg,coef,idx_ribosome_formation);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11) Dilution of mRNA, ribosome and enzymes.

% Dilution of RNA.
idx = cellfun(@(x) isequal({'RNA dilution'},x),model.subSystems,...
              'UniformOutput',false);

rna_rxns = model.rxns(cell2mat(idx));
rna_dil_rxns = rna_rxns(~contains(rna_rxns,'_generic_tRNA_'));

coef = mu/(k.deg_mrna + mu);

for i = 1:length(rna_dil_rxns)
    rxn_id = rna_dil_rxns{i};
    idx_dil = find(strcmp(model.rxns,rxn_id));
    rna_name = strrep(rxn_id,'R_dilution_mRNA_','');
    if contains(rxn_id,'sub')
        rna_met_name = strcat(rna_name,'_c');
        met_idx = strcmp(model.mets,rna_met_name);
        idx_cleavage = find(full(model.S(met_idx,:)) > 0);        

        if length(idx_cleavage) == 1
            fprintf(fptr,'CRNAdil%d: X%d - %.15f X%d = 0\n',...
                         i,idx_dil,coef,idx_cleavage);
        else
            for j = 1:length(idx_cleavage)
                idx = idx_cleavage(j);
            if j == 1
                eq = sprintf('X%d - %.15f X%d',idx_dil,coef,idx);
            else
                eq = sprintf('%s - %.15f X%d',eq,coef,idx);
            end
            end
            fprintf(fptr,'CRNAdil%d: %s = 0\n',i,eq);
        end
    else
        transcription_id = strcat('R_',rna_name);
        idx_transcription = find(strcmp(model.rxns,transcription_id));        
        fprintf(fptr,'CRNAdil%d: X%d - %.15f X%d = 0\n',...
                     i,idx_dil,coef,idx_transcription);
    end
end

% Dilution of enzymes.
idx = cellfun(@(x) isequal({'Enzyme dilution'},x),model.subSystems,...
              'UniformOutput',false);
enzyme_dil_rxns = model.rxns(cell2mat(idx));
enzyme_dil_rxns = enzyme_dil_rxns(contains(enzyme_dil_rxns,'_Enzyme'));

coef = mu/(k.deg_enzyme + mu);

for i = 1:length(enzyme_dil_rxns)
    rxn_id = enzyme_dil_rxns{i};
    idx_dil = find(strcmp(model.rxns,rxn_id));
    enzyme_name = strrep(rxn_id,'R_dilution_','');
    enzyme_formation_id = strcat('R_',enzyme_name);
    idx_formation = find(strcmp(model.rxns,enzyme_formation_id));
    fprintf(fptr,'Cenzymedil%d: X%d - %.15f X%d = 0\n',...
                 i,idx_dil,coef,idx_formation);
end
% Dilution of unmodeled protein.
ump_dil_idx = find(strcmp(model.rxns,'R_dilute_unmodeled_protein_to_biomass'));
ump_syn_idx = find(strcmp(model.rxns,'R_unmodeled_protein'));
coef = mu/(k.deg_enzyme + mu);
fprintf(fptr,'Cumpdil: X%d - %.15f X%d = 0\n',...
              ump_dil_idx,coef,ump_syn_idx);

% Dilution of 70S ribosome.

coef = mu/(k.deg_ribosome + mu);

ribosome_dil_rxn = 'R_dilution_ribosome_70S';
idx_dil = find(strcmp(model.rxns,ribosome_dil_rxn));
ribosome_formation_rxn = 'R_ribosome_70S';
idx_ribosome_formation = find(strcmp(model.rxns,ribosome_formation_rxn));
fprintf(fptr,'Cribodil: X%d - %.15f X%d = 0\n',...
             idx_dil,coef,idx_ribosome_formation);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12) Constraint of total protein and RNA.

dil_rxns = model.rxns(contains(model.rxns,'R_dilution'));

% %Sum of total protein and RNA.
% for i = 1:length(dil_rxns)
%     rxn_id = dil_rxns{i};
%     comp_name = strrep(rxn_id,'R_dilution_','');
%     idx = find(strcmp(model.rxns,rxn_id));
%     
%     if mod(i,200) == 0
%         sep = newline;
% 	else
%         sep = '';
%     end
%     
%     if strcmp(comp_name,'ribosome_70S')
%         MW = Info_ribosome.MW_protein + Info_ribosome.MW_RNA;
%         coeff = MW/1000;
%         if i == 1
%             eq = sprintf('%.15f X%d',coeff,idx);
%         else
%             eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
%         end
%         
%     elseif contains(comp_name,'_Enzyme')
%         MW = Info_enzyme.MW(contains(Info_enzyme.ID,comp_name));
%         coeff = MW/1000;
%         if i == 1
%             eq = sprintf('%.15f X%d',coeff,idx);
%         else
%             eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
%         end
%     
% 	elseif contains(comp_name,'generic_tRNA_')
%         trna_name = strrep(comp_name,'generic_','');
%         ge_trna_list = model.rxns(contains(model.rxns,'R_Generic_RNA_'));
%         trna_list = ge_trna_list(contains(ge_trna_list,trna_name));
%         trna_list = cellfun(@(x) strrep(x,'R_Generic_RNA_',''),...
%                                  trna_list,'UniformOutput',false);
%         trna_list = cellfun(@(x) x(1:min(strfind(x,trna_name))-2),...
%                                  trna_list,'UniformOutput',false);
%         MW_list = Info_tRNA.MW(ismember(Info_tRNA.ID,trna_list));
%         MW = mean(MW_list);
%         coeff = MW/1000;
%         if i == 1
%             eq = sprintf('%.15f X%d',coeff,idx);
%         else
%             eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
%         end
%         
%     elseif contains(comp_name,'mRNA_TU')
%         tu_name = strrep(comp_name,'mRNA_','');
%         MW = Info_mRNA.MW(strcmp(Info_mRNA.ID,tu_name));
%         coeff = MW/1000;       
%         if i == 1
%             eq = sprintf('%.15f X%d',coeff,idx);
%         else
%             eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
%         end
%     end
% end
% 
% fprintf(fptr,'CprotRNA: %s = %.15f\n',eq,mu*f);
% 
% %Total RNA constraint
% dil_rxns = dil_rxns(~contains(dil_rxns,'_Enzyme'));
% 
% for i = 1:length(dil_rxns)
%     rxn_id = dil_rxns{i};
%     comp_name = strrep(rxn_id,'R_dilution_','');
%     idx = find(strcmp(model.rxns,rxn_id));
%     
%     if mod(i,200) == 0
%         sep = newline;
% 	else
%         sep = '';
%     end
%     
%     if strcmp(comp_name,'ribosome_70S')
%         MW = Info_ribosome.MW_RNA;
%         coeff = MW/1000;
%         if i == 1
%             eq = sprintf('%.15f X%d',coeff,idx);
%         else
%             eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
%         end
%     
% 	elseif contains(comp_name,'generic_tRNA_')
%         trna_name = strrep(comp_name,'generic_','');
%         ge_trna_list = model.rxns(contains(model.rxns,'R_Generic_RNA_'));
%         trna_list = ge_trna_list(contains(ge_trna_list,trna_name));
%         trna_list = cellfun(@(x) strrep(x,'R_Generic_RNA_',''),...
%                                  trna_list,'UniformOutput',false);
%         trna_list = cellfun(@(x) x(1:min(strfind(x,trna_name))-2),...
%                                  trna_list,'UniformOutput',false);
%         MW_list = Info_tRNA.MW(ismember(Info_tRNA.ID,trna_list));
%         MW = mean(MW_list);
%         coeff = MW/1000;
%         if i == 1
%             eq = sprintf('%.15f X%d',coeff,idx);
%         else
%             eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
%         end
%         
%     elseif contains(comp_name,'mRNA_TU')
%         tu_name = strrep(comp_name,'mRNA_','');
%         MW = Info_mRNA.MW(strcmp(Info_mRNA.ID,tu_name));
%         coeff = MW/1000;       
%         if i == 1
%             eq = sprintf('%.15f X%d',coeff,idx);
%         else
%             eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
%         end
%     end
% end
% 
% fprintf(fptr,'CtotRNA: %s = %.15f\n',eq,mu*0.56*(mu+1.25)/(mu+12.51));

%Total protein constraint
for i = 1:length(dil_rxns)
    rxn_id = dil_rxns{i};
    comp_name = strrep(rxn_id,'R_dilution_','');
    idx = find(strcmp(model.rxns,rxn_id));
    
    if mod(i,200) == 0
        sep = newline;
	else
        sep = '';
    end
    
    if strcmp(comp_name,'ribosome_70S')
        MW = Info_ribosome.MW_protein;
        coeff = MW/1000;
        if i == 1
            eq = sprintf('%.15f X%d',coeff,idx);
        else
            eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
        end
    
    elseif contains(comp_name,'_Enzyme')
        MW = Info_enzyme.MW(contains(Info_enzyme.ID,comp_name));
        coeff = MW/1000;
        if i == 1
            eq = sprintf('%.15f X%d',coeff,idx);
        else
            eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
        end
    end
end

fprintf(fptr,'Ctotprot: %s = %.15f\n',eq,mu*f);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 13) Constraint of membrane proteins (transporters).

% Transporters = {};
% idx = cellfun(@(x) isequal({'Metabolism'},x),model.subSystems,'UniformOutput',false);
% M_rxns_idx = find(cell2mat(idx));
% for i = 1:length(M_rxns_idx)
%     rxnidx = M_rxns_idx(i);
% 	metlist = model.mets(model.S(:,rxnidx) ~= 0);
% 	if find(cell2mat(cellfun(@(x) strcmp(x(length(x)-1:end),'_e'),metlist,'UniformOutput',false)),1) > 0
%         enzyme_id = model.grRules{rxnidx};
%         if ~isempty(enzyme_id)
%             Transporters = [Transporters;{enzyme_id}];
%         end
% 	end
% end

Transporters = {'M_GLCpts_1_Enzyme_c';'M_GLCpts_2_Enzyme_c';'M_GLCpermease_fwd_Enzyme_c'};

for i = 1:length(Transporters)
	comp_name = Transporters{i};
	rxn_id = strcat('R_dilution_',comp_name(1:length(comp_name)-2));
	idx = find(strcmp(model.rxns,rxn_id));

	if mod(i,200) == 0
        sep = newline;
    else
        sep = '';
	end

	MW = Info_enzyme.MW(contains(Info_enzyme.ID,comp_name));
	coeff = MW/1000;
	if i == 1
        eq = sprintf('%.15f X%d',coeff,idx);
    else
        eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
	end
end

fprintf(fptr,'Ctransporter: %s <= %.15f\n',eq,mu*0.46*f_transporter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set lower and upper bounds.

fprintf(fptr,'Bounds\n');

for i = 1:length(model.rxns)
	if model.ub(i) >= 100
        fprintf(fptr,'%f <= X%d <= +infinity\n',model.lb(i),i);
    else
        fprintf(fptr,'%f <= X%d <= %f\n',model.lb(i),i,model.ub(i));
	end
end

fprintf(fptr,'End\n');
fclose(fptr);
