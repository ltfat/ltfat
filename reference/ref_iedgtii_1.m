function f=ref_iedgtii_1(c,g,a,M)
%REF_EDGTII_1   Reference Inverse Even DGT type II by DGT
%   Usage  c=ref_edgt(f,g,a,M);
%
%   The input window must be odd-centered of length 2L.
%
%   a must be divisable by 2.

L=size(g,1)/2;
W=size(c,2);

N=L/a;

clong=zeros(M,2*N,W);

cr=reshape(c,M,N,W);
% Copy the first half unchanged
clong(:,1:N,:)=cr;

% Copy the non modulated coefficients.
clong(1,N+1:2*N,:)=cr(1,N:-1:1,:);

% Copy the modulated coefficients.
clong(2:M,N+1:2*N,:)=-cr(M:-1:2,N:-1:1,:);

clong=reshape(clong,2*M*N,W);

fdouble=ref_igdgt(clong,g,a,M,.5,0,floor(a/2));

f=fdouble(1:L,:);


