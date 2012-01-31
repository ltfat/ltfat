function c=ref_dft(f)
%REF_DFT  Reference Discrete Fourier Transform
%   Usage:  c=ref_dft(f);
%
%   REF_DFT(f) computes a normalized discrete Fourier transform of 
%   f, so this is not the same as FFT(f).

L=size(f,1);
W=size(f,2);

% Create weights.
w=sqrt(1/L);

% Create transform matrix.
F=zeros(L);

for m=0:L-1
  for n=0:L-1
    F(m+1,n+1)=w*exp(2*pi*i*m*n/L);
  end;
end;

% Compute coefficients.
c=F'*f;


