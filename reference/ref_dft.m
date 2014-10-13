function c=ref_dft(f)
%REF_DFT  Reference Discrete Fourier Transform
%   Usage:  c=ref_dft(f);
%
%   REF_DFT(f) computes the unitary discrete Fourier transform of f.
%
%   AUTHOR: Jordy van Velthoven

L = length(f);
c = zeros(L,1);

for k=0:L-1
  for l=0:L-1
    c(k+1) = c(k+1) + f(l+1) * exp(-2*pi*i*k*l/L);
  end;
end;

c = c/sqrt(L);


