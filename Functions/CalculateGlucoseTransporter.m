%% CalculateGlucoseTransporter
function value = CalculateGlucoseTransporter(model,sol_full,Info_enzyme)
% The unit is g/gCDW.
enzymes = {'M_GLCpts_1_Enzyme_c';
           'M_GLCpts_2_Enzyme_c';
           'M_GLCt2_fwd_Enzyme_c';
           'M_GLCt2_rvs_Enzyme_c'};
[~, b] = size(sol_full);
value = zeros(1,b);
for i = 1:b
    flux_tmp = sol_full(:,i);
    value(1,i) = sum(CalculateEnzymeWeight(model,enzymes,flux_tmp,Info_enzyme));
end
