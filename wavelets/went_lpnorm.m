function E = went_lpnorm(x,p)
%WENT_LPNORM Lp-norm non-additive entropy
%   Usage:  E = went_lpnorm(x,p)
%
%   `E = went_lpnorm(x,p)` calculates entropy of the input *x* using formula
%   XXX
%

assert(0<p&&p<2,'%s: Parameter p have to be in the range ]0,2[',upper(mfilename));
%u = x./norm(x);

%E = sum(abs(x(:)).^p)^(1/p);
E = sum(abs(x(:)).^p);