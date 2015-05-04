function R=ref_drihaczekdist(f)
%REF_DRIHACZEKDIST  Reference discrete Rihaczek distribution
%   Usage:  R=ref_drihaczekdist(f)
%
%   REF_DRIHACZEKDIST(f) computes the discrete Rihaczek distribution of f.

% AUTHOR: Jordy van Velthoven

L = length(f);
c = fft(f);


for m = 0:L-1
  for n=0:L-1
    R(m+1, n+1) = (f(n+1).' * conj(c(m+1))) .* exp(-2*pi*i*m*n/L);
   end
end
