function c=ref_dftiv(f)
%REF_DFT  Reference Discrete Fourier Transform Type IV
%   Usage:  c=ref_dftiv(f);
%
%   This is highly experimental!

L=size(f,1);
W=size(f,2);

% Create weights.
w=sqrt(1/L);

% Create transform matrix.
F=zeros(L);

for m=0:L-1
  for n=0:L-1
    F(m+1,n+1)=w*exp(2*pi*i*(m+.5)*(n+.5)/L);
  end;
end;

% Compute coefficients.
c=F'*f;


