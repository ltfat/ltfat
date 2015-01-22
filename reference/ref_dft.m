function c=ref_dft(f)
%REF_DFT  Reference Discrete Fourier Transform
%   Usage:  c=ref_dft(f);
%
%   REF_DFT(f) computes the unitary discrete Fourier transform of f.
%
%   AUTHOR: Jordy van Velthoven

L=size(f,1);
W=size(f,2);
c= zeros(L,W);

for w=1:W
for k=0:L-1
  for l=0:L-1
    c(k+1,w) = c(k+1,w) + f(l+1,w) * exp(-2*pi*i*k*l/L);
  end;
end;
end

c = c./sqrt(L);


