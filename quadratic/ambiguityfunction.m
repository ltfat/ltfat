function A = ambiguityfunction(f,g)
%AMBIGUITYFUNCTION Ambiguity function
%   Usage: A = ambiguityfunction(f);
%          A = ambiguityfunction(f,g);
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

upfname = upper(mfilename);
complainif_notenoughargs(nargin, 1, upfname);
complainif_toomanyargs(nargin, 2, upfname);

[f,~,W]=comp_sigreshape_pre(f,upfname);

if W>1
   error('%s: Only one-dimensional vectors can be processed.',upfname); 
end

if (nargin == 1)
  if isreal(f)
    z1 = comp_fftanalytic(f);
  else
    z1 = f;
  end
  z2 = z1;

elseif (nargin == 2)
  [g,~,W]=comp_sigreshape_pre(g,upfname);
  if W>1
     error('%s: Only one-dimensional vectors can be processed.',upfname); 
  end

  if ~all(size(f)==size(g))
  	error('%s: f and g must have the same length.', upfname);
  end;

  if xor(isreal(f), isreal(g))
      error('%s: One input is real, the other one must be real too. ',...
            upfname);
  end

  if isreal(f) || isreal(g)
    z1 = comp_fftanalytic(f);
    z2 = comp_fftanalytic(g);
  else
    z1 = f;
    z2 = g;
  end;
end

R = comp_instcorrmat(z1, z2);

A = fftshift(fft2(fft(R)));


