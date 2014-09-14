function r=lxcorr(f,g,varargin)
%LXCORR  Linear cross-correlation
%   Usage:  h=lxcorr(f,g)
%
%   `lxcorr(f)` computes the linear cross-correlation of the input
%   signal *f*. The linear cross-correlation is defined by
%
%   ..          L-1
%      h(l+1) = sum f(k+1) * conj(g(k-l+1))
%               k=0
%
%   .. math:: h\left(l+1\right)=\sum_{k=0}^{L-1}f\left(k+1\right)\overline{g\left(k-l+1\right)}

%
%   `lxcorr(f,'normalize')` does the same, but normalizes the output by
%   the product of the $l^2$-norm of *f* and *g*.
%
%   See also: paucorr, lxcorr

%   AUTHOR: Jordy van Velthoven

definput.flags.type={'nonormalize','normalize'};

flags = ltfatarghelper({},definput,varargin);

h = lconv(f, g, 'r');

if flags.do_normalize
  h = h/(norm(f)*norm(g));  
end
