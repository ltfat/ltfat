function c=ref_dctii(f)
%REF_DCTII  Reference Discrete Consine Transform type II
%   Usage:  c=ref_dctii(f);
%
%

L=size(f,1);
W=size(f,2);

% Create weights.
w=ones(L,1);
w(1)=1/sqrt(2);
w=w*sqrt(2/L);

% Create transform matrix.
F=zeros(L);

for m=0:L-1
  for n=0:L-1
    F(m+1,n+1)=w(n+1)*cos(pi*n*(m+.5)/L);
  end;
end;

% Compute coefficients.
c=F'*f;


