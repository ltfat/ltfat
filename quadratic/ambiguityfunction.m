function A = ambiguityfunction(f,g);
%AMBIGUITYFUNCTION Ambiguity function
%   Usage A = ambiguityfunction(f);
%         A = ambiguityfunction(f,g);
%
%   Input parameters:
%         f,g    : Input vector(s).
%
%   Output parameters:
%         A      : ambiguity function
%
%   `ambiguityfunction(f)` computes the (symmetric) ambiguity function of f.
%   The ambiguity function is computed as the two-dimensional Fourier transform
%   of the Wigner-Ville distribution |wignervilledist|.
%
%   **WARNING**: The quadratic time-frequency distributions are highly 
%   redundant. For an input vector of length L, the quadratic time-frequency
%   distribution will be a $L \times L$ matrix.

% AUTHOR: Jordy van Velthoven
% TESTING: TEST_AMBIGUITYFUNCTION
% REFERENCE: REF_AMBIGUITYFUNCTION

complainif_notenoughargs(nargin, 1, 'DAF');

if (nargin == 1)

  [f,~,Lf,W,~,permutedsize,order]=assert_sigreshape_pre(f,[],[],upper(mfilename));
  
  if isreal(f)
    z1 = comp_fftanalytic(f, Lf);
    z2 = z1;
  else
    z1 = f;
    z2 = z1;
  end
 
elseif (nargin == 2)
  [f,~,Lf,W,~,permutedsize,order]=assert_sigreshape_pre(f,[],[],upper(mfilename));
  [g,~,Lg,W,~,permutedsize,order]=assert_sigreshape_pre(g,[],[],upper(mfilename));

  if ~all(size(f)==size(g))
  	error('%s: f and g must have the same length.', upper(mfilename));
  end;
  
  if isreal(f) || isreal(g)
    z1 = comp_fftanalytic(f, Lf);
    z2 = comp_fftanalytic(g, Lg);
  else
    z1 = f;
    z2 = g;
  end;
end

R = comp_instcorrmat(z1, z2, Lf);

A = fftshift(fft2(fft(R)));