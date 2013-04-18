function g=comp_pchirp(L,n)
%COMP_PCHIRP  Compute periodic chirp
%   Usage:  g=comp_pchirp(L,n);
%
%   `pchirp(L,n)` returns a periodic, discrete chirp of length *L* that
%   revolves *n* times around the time-frequency plane in frequency. *n* must be
%   an integer number.
%
%   This is a computational routine. Do not call it unless you have
%   verified the input parameters.

%   AUTHOR : Peter L. Søndergaard
%   TESTING: OK
%   REFERENCE: OK

l= (0:L-1).';
X = mod(mod(mod(n*l,2*L).*l,2*L)*(L+1),2*L);
g = exp(pi*1i*X/L);
