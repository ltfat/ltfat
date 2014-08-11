function h=pxcorr(f,g,varargin)
%PXCORR  Periodic cross correlation
%   Usage:  h=pxcorr(f,g)
%
%   `pxcorr(f,g)` computes the periodic cross correlation of the input
%   signals *f* and *g*. The cross correlation is defined by
%
%   ..          L-1
%      h(l+1) = sum f(k+1) * conj(g(k-l+1))
%               k=0
%
%   .. math:: h\left(l+1\right)=\sum_{k=0}^{L-1}f\left(k+1\right)\overline{g\left(k-l+1\right)}
%
%   In the above formula, $k-l$ is computed modulo $L$.
%
%   `pxcorr(f,g,'normalize')` does the same, but normalizes the output by
%   the product of the $l^2$-norm of *f* and *g*.
%
%   See also: dft, pfilt, involute

%   AUTHOR: Peter L. Søndergaard, Jordy van Velthoven

definput.flags.type={'nonormalize','normalize'};

flags = ltfatarghelper({},definput,varargin);

h = pconv(f, g, 'r');

if flags.do_normalize
  h = h/(norm(f)*norm(g));  
end


