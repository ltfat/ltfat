function W = wignervilledist(f,varargin)
%WIGNERVILLEDIST Wigner-Ville distribution
%   Usage: W = wignervilledist(f);
%          W = wignervilledist(f, g);
%
%   Input parameters:
%         f,g      : Input vector(s)
%
%   Output parameters:
%         w      : Wigner-Ville distribution
%
%   `wignervilledist(f)` computes the Wigner-Ville distribution of the vector f. The
%   Wigner-Ville distribution is computed by
%
%   .. math:: W\left( n+1,k+1 \right)\; =\; \sum_{m=0}^{L-1}{R\left( n+1,m+1 \right)e^{-i2\pi mk/L}},
%
%   where $R(n,m)$ is the instantaneous correlation matrix given by
%
%   .. math:: R\left( n,m \right)\; =\; z\left( n+m \right)\overline{z\left( n-m \right)},
%
%   where $m \in {-L/2,\ldots, L/2 - 1}$, and where $z$ is the analytical representation of
%   $f$, when $f$ is real-valued.
%
%   `wignervilledist(f,g)` computes the cross-Wigner-Ville distribution of f and g.
%
%   **WARNING**: The quadratic time-frequency distributions are highly
%   redundant. For an input vector of length L, the quadratic time-frequency
%   distribution will be a $L \times L$ matrix.

% AUTHOR: Jordy van Velthoven
% TESTING: TEST_WIGNERVILLEDIST
% REFERENCE: REF_WIGNERVILLEDIST

upfname = upper(mfilename);
complainif_notenoughargs(nargin, 1, upfname);

definput.flags.complex = {'asinput','complex'};
definput.keyvals.g=[];
[flags,kv,g]=ltfatarghelper({'g'},definput,varargin);

[f,~,W]=comp_sigreshape_pre(f,upfname);
if W>1
   error('%s: Only one-dimensional vectors can be processed.',upfname); 
end

if isempty(g)
  if isreal(f) && ~flags.do_complex
    z1 = comp_fftanalytic(f);
  else
    z1 = f;
  end
  
  z2 = z1;
  R = comp_instcorrmat(z1, z2);

  W = real(fft(R));

else
  [g,~,W]=comp_sigreshape_pre(g,upfname);
  
  if W>1
    error('%s: Only one-dimensional vectors can be processed.',upfname); 
  end

  if ~all(size(f)==size(g))
  	error('%s: f and g must have the same length.', upper(mfilename));
  end;
  
  if xor(isreal(f), isreal(g))
      error('%s: One input is real, the other one must be real too. ',...
            upfname);
  end

  if isreal(f) || isreal(g) && ~flags.do_complex
    z1 = comp_fftanalytic(f);
    z2 = comp_fftanalytic(g);
  else
    z1 = f;
    z2 = g;
  end;

  R = comp_instcorrmat(z1, z2);

  W = fft(R);
end
