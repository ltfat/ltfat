function f=ref_idft(c)
%REF_IDFT  Reference Inverse Discrete Fourier Transform
%   Usage:  f=ref_dft(c);
%
%   REF_IDFT(f) computes the unitary discrete Fourier transform of the 
%   coefficient c.
%
%   AUTHOR: Jordy van Velthoven

L = length(c);
f = zeros(L,1);

for l=0:L-1
  for k=0:L-1
    f(l+1) = f(l+1) + c(k+1) * exp(2*pi*i*k*l/L);
  end;
end;

f = f/sqrt(L);
