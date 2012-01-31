function f=ref_irdgt_1(c,g,a,M)
%REF_IRDGT_1  Reference Inverse Real DGT by fac. and IRDFT
%   Usage:  c=ref_rdgt_1(f,g,a,M);
%
%   Compute the factorization and use IRDFT

L=size(g,1);
R=size(g,2);
W=size(c,2);

%M=L/b;
N=L/a;

% Apply ifft to the coefficients.
c=ref_irdft(reshape(c,M,N*R*W));

gf = comp_wfac(g,a,M);      
f = comp_idgt_fw(c,gf,L,a,M);


