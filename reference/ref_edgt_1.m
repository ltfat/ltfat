function c=ref_edgt_1(f,g,a,M)
%REF_EDGT_1   Reference Even Discrete Gabor transform by DGT
%   Usage  c=ref_edgt(f,g,a,M);
%
%   The input window must be odd-centered of length 2L.
%
%   M must be even

L=size(f,1);
W=size(f,2);

N=L/a;
M=L/b;

clong=ref_dgtii([f;flipud(f)],g,a,2*b);

c=zeros(M*N,W);

% The coefficient array is stacked from:
% - The first M/2 coefficients of the first time shift.
% - The body of the coefficients 
% - The first M/2 coefficients of the last time shift.

c(:,:)=[clong(1:M/2,:); ...
        clong(M+1:M*N+M/2,:)];

% Scale the first coefficients correctly
c(1,:)=c(1,:)/sqrt(2);
c(M*(N-1)+M/2+1,:)=c(M*(N-1)+M/2+1,:)/sqrt(2);


