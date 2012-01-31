function c=ref_dctiv(f)
%REF_DCTIV  Reference Discrete Consine Transform type IV
%   Usage:  c=ref_dctiv(f);
%

L=size(f,1);
W=size(f,2);

% Create weights.
w=sqrt(2/L);

% Create transform matrix.
F=zeros(L);

for m=0:L-1
  for n=0:L-1
    F(m+1,n+1)=w*cos(pi*(n+.5)*(m+.5)/L);
  end;
end;

% Compute coefficients.
c=F*f;


