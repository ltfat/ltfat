function qgc = gc_mu_uquant(gc,Nbits,mu)
%
% Uniform quantization of Gabor coefficients after mu-law companding
%   qgc = gc_mu_uquant(gc,Nbits,mu)
%

gc_r = real(gc);
gc_i = imag(gc);

[mu_gc_r,weight_r] = mulawencode(gc_r,mu);
[mu_gc_i,weight_i] = mulawencode(gc_i,mu);
qgc_r = uquant(mu_gc_r,Nbits);
qgc_i = uquant(mu_gc_i,Nbits);
qgc_r = mulawdecode(qgc_r,mu,weight_r);
qgc_i = mulawdecode(qgc_i,mu,weight_i);
qgc = qgc_r + i* qgc_i;


