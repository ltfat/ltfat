function c=ref_rdgt_1(f,g,a,M)
%REF_RDGT_1  Reference Real DGT by fac. and RDFT
%   Usage:  c=ref_rdgt_1(f,g,a,M);
%
%   Compute the factorization and use RDFT

L=size(f,1);
W=size(f,2);
R=size(g,2);

N=L/a;

gf=comp_wfac(g,a,M);

% Compute the window application and the DFT modulation.
c=zeros(M*N*R,W);
c(:)=ref_rdft(comp_dgt_fw(f,gf,L,a,M));


