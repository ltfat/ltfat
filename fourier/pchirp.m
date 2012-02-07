function g=pchirp(L,n)
%PCHIRP  Periodic chirp
%   Usage:  g=pchirp(L,n);
%
%   `pchirp(L,n)` returns a periodic, discrete chirp of length *L* that
%   revolves *n* times around the time-frequency plane in frequency. *n* must be
%   a whole number.
%
%   To get a chirp that revolves around the time-frequency plane in time,
%   use ::
%
%     dft(pchirp(L,N));  
%
%   The chirp is computed by:
%   
%   ..  g(l+1) = exp(pi*i*n*l^2/L) for l=0,...,L-1
%
%   .. math:: g\left(l+1\right)=e^{\pi inl^{2}/L},\quad l=0,\ldots,L-1
%
%   The chirp has absolute value 1 everywhere. To get a chirp with unit
%   $l^2$-norm, divide the chirp by $\sqrt(L)$.
%
%   See also: dft, expwave
%
%   References: fehakamane06

%   AUTHOR : Peter Soendergaard
%   TESTING: OK
%   REFERENCE: OK

error(nargchk(2,2,nargin));

% Compute normalized chirp
g=(exp((0:L-1).^2/L*pi*i*n)/sqrt(L)).';

