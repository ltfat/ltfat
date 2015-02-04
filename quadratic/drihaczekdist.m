function r = drihaczekdist(f);
%DRIHACZEKDIST discrete Rihaczek distribution
%   Usage r = drihaczekdist(f);
%
%
%   `drihaczekdist(f)` computes a discrete Rihaczek distribution.
%   The discrete Rihaczek distribution is computed by
%
%   .. math:: r\left( k+1,\; l+1 \right)\; =\; f\left( l+1 \right)\; \overline{c\left( k+1 \right)}e^{-2\pi ikl/L}
%
%   where $k, l=0,\ldots,L-1$ and $c$ is the Fourier transform of $f$.
%
%   **WARNING**: The quadratic time-frequency distributions are highly 
%   redundant. For an input vector of length L, the quadratic time-frequency
%   distribution will be a $L \times L$ matrix.
%
%   AUTHOR: Jordy van Velthoven
%   TESTING: TEST_DRIHACZEKDIST
%   REFERENCE: REF_DRIHACZEKDIST

complainif_notenoughargs(nargin, 1, 'DRIHACZEKDIST');

[f,~,Ls,W,~,permutedsize,order]=assert_sigreshape_pre(f,[],[],upper(mfilename));

c = dgt(f, f, 1, Ls);

r = dsft(c);
