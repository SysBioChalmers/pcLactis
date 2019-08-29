%% SetParameters 
function k = SetParameters(mu)

k = struct();

% Kinetic parameters
k.deg_mrna = 6.25*mu+1.54; %mRNA degradation constant (/h)
k.deg_enzyme = 0.01*mu; %protein degradation constant (/h)
k.deg_ribosome = 0.01*mu; %ribosome degradation constant (/h)
% Catalytic rates.
k.cplx_rnap = 176.46*mu/(1.147+mu); %nucleotide molecules/complex molecule/s
k.cplx_mrnadeg = 100*1635*1.43*mu; %nucleotide molecules/complex molecule/s
k.trna = 6.65*mu/(1.147+mu); %amino acid molecules/tRNA molecule/s
k.mrna = 2.04*mu/(1.147+mu); %protein molecules/TU molecule/s
k.ribo = 58.82*mu/(1.147+mu); %amino acid molecules/ribosome molecule/s
k.protease = 100*351*1.43*mu; %amino acid molecules/complex molecule/s

% Convert to 1/h
k.cplx_rnap = 3600*k.cplx_rnap; %nucleotide molecules/complex molecule/h
k.cplx_mrnadeg = 3600*k.cplx_mrnadeg; %nucleotide molecules/complex molecule/h
k.trna = 3600*k.trna; %amino acid molecules/tRNA molecule/h
k.mrna = 3600*k.mrna; %protein molecules/TU molecule/h
k.ribo = 3600*k.ribo; %amino acid molecules/ribosome molecule/h
k.protease = 3600*k.protease; %amino acid molecules/complex molecule/h