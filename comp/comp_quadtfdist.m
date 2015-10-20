function p = comp_quadtfdist(f, q);;
% Comp_QUADTFDIST Compute quadratic time-frequency distribution
%   Usage p = comp_quadtfdist(f, q);;
%
%   Input parameters:
%         f      : Input vector
%	  q	 : Kernel
%
%   Output parameters:
%         p      : Quadratic time-frequency distribution
%

% AUTHOR: Jordy van Velthoven

if isreal(f)
 z = comp_fftanalytic(f);
else
 z = f;
end;

R = comp_instcorrmat(z,z);

c = ifft2(fft2(R).*fft2(q));

p = fft(c);


