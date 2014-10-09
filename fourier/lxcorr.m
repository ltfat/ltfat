function h=lxcorr(f,g,varargin)
%LXCORR  Linear crosscorrelation
%   Usage:  h=lxcorr(f,g)
%
%   `lxcorr(f)` computes the linear crosscorrelation of the input signal *f* and *g*. 
%   The linear cross-correlation is computed by
%
%   ..          Lh-1
%      h(l+1) = sum f(k+1) * conj(g(k-l+1))
%               k=0
%
%   .. math:: h\left(l+1\right)=\sum_{k=0}^{L-1}f\left(k+1\right)\overline{g\left(k-l+1\right)}
%
%   with $L_{h} = L_{f} + L_{g} - 1$ where $L_{f}$ and $L_{g}$ are the lengths of *f* and *g*, 
%   respectively.
%
%   `lxcorr(f,'normalize')` does the same, but normalizes the output by
%   the product of the $l^2$-norm of *f* and *g*.
%
%   See also: pxcorr, lconv

%   AUTHOR: Jordy van Velthoven

definput.flags.type={'nonormalize','normalize'};

flags = ltfatarghelper({},definput,varargin);

h = lconv(f, g, 'r');

if flags.do_normalize
  h = h/(norm(f)*norm(g));  
end
