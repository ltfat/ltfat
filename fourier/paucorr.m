function r=paucorr(f,varargin)
%PAUCORR  Periodic autocorrelation
%   Usage:  h=paucorr(f)
%
%   `paucorr(f)` computes the periodic autocorrelation of the input
%   signal *f*. The autocorrelation is defined by
%
%   ..          L-1
%      r(l+1) = sum f(k+1) * conj(f(k-l+1))
%               k=0
%
%   .. math:: r\left(l+1\right)=\sum_{k=0}^{L-1}f\left(k+1\right)\overline{f\left(k-l+1\right)}
%
%   In the above formula, $k-l$ is computed modulo $L$.
%
%   `paucorr(f,'normalize')` does the same, but normalizes the $l^2$-norm of *f* to be 1.
%
%   See also: pxcorr, pconv, dft

%   AUTHOR: Jordy van Velthoven

r = pxcorr(f,f,varargin{:});
