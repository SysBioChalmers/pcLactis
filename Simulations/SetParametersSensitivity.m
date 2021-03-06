%% SetParameters 
function k = SetParametersSensitivity(mu,ID,fold)

k = struct();

% Kinetic parameters
k.deg_mrna = 6.25*mu+1.54; %mRNA degradation constant (/h)
k.deg_enzyme = 0.01*mu; %protein degradation constant (/h)
k.deg_ribosome = 0.01*mu; %ribosome degradation constant (/h)

% Catalytic rates.

k.cplx_rnap = 116.4*mu/(0.55+mu); %nucleotide molecules/complex molecule/s
k.trna = 4.1*mu/(0.55+mu); %amino acid molecules/tRNA molecule/s
k.mrna = 1.45*mu/(0.55+mu); %protein molecules/TU molecule/s
k.ribo = 38.8*mu/(0.55+mu); %amino acid molecules/ribosome molecule/s

k.cplx_mrnadeg = 100*1591; %nucleotide molecules/complex molecule/s
k.protease = 100*350; %amino acid molecules/complex molecule/s

% Convert to 1/h
k.cplx_rnap = 3600*k.cplx_rnap; %nucleotide molecules/complex molecule/h
k.cplx_mrnadeg = 3600*k.cplx_mrnadeg; %nucleotide molecules/complex molecule/h
k.trna = 3600*k.trna; %amino acid molecules/tRNA molecule/h
k.mrna = 3600*k.mrna; %protein molecules/TU molecule/h
k.ribo = 3600*k.ribo; %amino acid molecules/ribosome molecule/h
k.protease = 3600*k.protease; %amino acid molecules/complex molecule/h

if ismember(ID,'mRNA_degradation_constant')
    k.deg_mrna = k.deg_mrna * fold;
elseif ismember(ID,'Protein_degradation_constant')
    k.deg_enzyme = k.deg_enzyme * fold;
elseif ismember(ID,'Ribosome_degradation_constant')
    k.deg_ribosome = k.deg_ribosome * fold;
elseif ismember(ID,'RNA_polymerase_catalytic_rate')
    k.cplx_rnap = k.cplx_rnap * fold;
elseif ismember(ID,'mRNA_degradation_complex_catalytic_rate')
    k.cplx_mrnadeg = k.cplx_mrnadeg * fold;
elseif ismember(ID,'tRNA_catalytic_rate')
    k.trna = k.trna * fold;
elseif ismember(ID,'mRNA_catalytic_rate')
    k.mrna = k.mrna * fold;
elseif ismember(ID,'Ribosomal_catalytic_rate')
    k.ribo = k.ribo * fold;
elseif ismember(ID,'Protease_complex_catalytic_rate')
    k.protease = k.protease * fold;
end

