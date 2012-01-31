function c=ref_edgtii_1(f,g,a,M)
%REF_EDGTII_1   Reference Even Discrete Gabor transform type II by DGT
%   Usage  c=ref_edgt(f,g,a,M);
%
%   If a is even, then the input window must be odd-centered of length 2L.
%   
%   If a is odd, then the input window must be even-centered of length 2L.
%

L=size(f,1);
W=size(f,2);

N=L/a;

clong=ref_gdgt([f;conj(flipud(f))],g,a,M,.5,0,floor(a/2));

c=clong(1:M*N,:);


