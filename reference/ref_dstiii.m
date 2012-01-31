function c=ref_dstiii(f)
%REF_DSTIII  Reference Discrete Sine Transform type III
%   Usage:  c=ref_dstiii(f);
%
%   This is the inverse of DSTII

L=size(f,1);
W=size(f,2);

% Create weights.
w=ones(L,1);
w(L)=1/sqrt(2);
w=w*sqrt(2/L);

% Create transform matrix.
F=zeros(L);

for m=0:L-1
  for n=0:L-1
    F(m+1,n+1)=w(n+1)*sin(pi*(n+1)*(m+.5)/L);
  end;
end;

% Compute coefficients.
c=F*f;


