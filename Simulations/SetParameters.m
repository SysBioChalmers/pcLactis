%% SetParameters 
function k = SetParameters(mu)

k = struct();

% Kinetic parameters
k.deg_mrna = 6.25*mu+1.54; %mRNA degradation constant (/h)
k.deg_enzyme = 0.1*mu; %protein degradation constant (/h)
k.deg_ribosome = 0.1*mu; %ribosome degradation constant (/h)

% Catalytic rates.
% k.cplx_rnap = 176.46*mu/(1.147+mu); %nucleotide molecules/complex molecule/s
% k.trna = 6.65*mu/(1.147+mu); %amino acid molecules/tRNA molecule/s
% k.mrna = 2.04*mu/(1.147+mu); %protein molecules/TU molecule/s
% k.ribo = 58.82*mu/(1.147+mu); %amino acid molecules/ribosome molecule/s

% k.cplx_rnap = 108.3*mu/(0.84+mu); %nucleotide molecules/complex molecule/s
% k.trna = 4.08*mu/(0.84+mu); %amino acid molecules/tRNA molecule/s
% k.mrna = 1.246*mu/(0.84+mu); %protein molecules/TU molecule/s
% k.ribo = 36.1*mu/(0.84+mu); %amino acid molecules/ribosome molecule/s

% k.cplx_rnap = 99*mu/(0.55+mu); %nucleotide molecules/complex molecule/s
% k.trna = 3.6*mu/(0.55+mu); %amino acid molecules/tRNA molecule/s
% k.mrna = 1.4*mu/(0.55+mu); %protein molecules/TU molecule/s
% k.ribo = 33*mu/(0.55+mu); %amino acid molecules/ribosome molecule/s

% k.cplx_rnap = 117.6*mu/(0.55+mu); %nucleotide molecules/complex molecule/s
% k.trna = 3.77*mu/(0.55+mu); %amino acid molecules/tRNA molecule/s
% k.mrna = 1.45*mu/(0.55+mu); %protein molecules/TU molecule/s
% k.ribo = 39.2*mu/(0.55+mu); %amino acid molecules/ribosome molecule/s

k.cplx_rnap = 116.4*mu/(0.55+mu); %nucleotide molecules/complex molecule/s
k.trna = 4.1*mu/(0.55+mu); %amino acid molecules/tRNA molecule/s
k.mrna = 1.45*mu/(0.55+mu); %protein molecules/TU molecule/s
k.ribo = 38.8*mu/(0.55+mu); %amino acid molecules/ribosome molecule/s

k.cplx_mrnadeg = 100*1591*1.43*mu; %nucleotide molecules/complex molecule/s
k.protease = 100*350*1.43*mu; %amino acid molecules/complex molecule/s

% Convert to 1/h
k.cplx_rnap = 3600*k.cplx_rnap; %nucleotide molecules/complex molecule/h
k.cplx_mrnadeg = 3600*k.cplx_mrnadeg; %nucleotide molecules/complex molecule/h
k.trna = 3600*k.trna; %amino acid molecules/tRNA molecule/h
k.mrna = 3600*k.mrna; %protein molecules/TU molecule/h
k.ribo = 3600*k.ribo; %amino acid molecules/ribosome molecule/h
k.protease = 3600*k.protease; %amino acid molecules/complex molecule/h