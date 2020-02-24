%% CalculateProteinCost
%   This function will calculate protein cost for each reaction and the
%   total cost for the pathway
function [protcost,totprotcost] = CalculateProteinCost(rxnslist,fluxlist,Info_enzyme,kcat_glc,enzymelist,kcatlist)

protcost = zeros(length(rxnslist),1);
for i = 1:length(rxnslist)
    rxnid = rxnslist{i};
    ezmid = strcat(rxnid(3:end),'_');
    flux = fluxlist(i);
    if ismember({rxnid},{'R_M_GLCpts'})
        kcat_tmp = kcat_glc * 3600;%/h
        mw_tmp = Info_enzyme.MW(ismember(Info_enzyme.ID,'M_GLCpts_2_Enzyme_c'));
        protcost_tmp = mw_tmp/kcat_tmp;
    else
        if any(contains(enzymelist,ezmid))
            enzymes_tmp = enzymelist(contains(enzymelist,ezmid));
            kcats_tmp = kcatlist(contains(enzymelist,ezmid));
            if flux < 0
                idx_tmp = contains(enzymes_tmp,'_rvs');
            else
                idx_tmp = ~contains(enzymes_tmp,'_rvs');
            end
            enzymes_tmptmp = enzymes_tmp(idx_tmp);
            kcats_tmptmp = kcats_tmp(idx_tmp);
            [~, b] = ismember(enzymes_tmptmp,Info_enzyme.ID);
            mws_tmp = Info_enzyme.MW(b);
            pcs_tmp = mws_tmp./kcats_tmptmp;
            protcost_tmp = min(pcs_tmp);
        else
            protcost_tmp = 0;
        end
    end
    protcost(i) = protcost_tmp; %gprotein/gCDW per flux(mol/gCDW/h)
end
totprotcost = sum(protcost.*abs(fluxlist))/1000; %gprotein/gCDW per flux(mmol/gCDW/h) of glucose


