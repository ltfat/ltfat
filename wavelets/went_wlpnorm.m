function E = went_wlpnorm(x,p)
%WENT_WLPNORM Weak Lp-norm non-additive entropy
%   Usage:  E = went_wlpnorm(x,p)
%
%   `E = went_wlpnorm(x,p)` calculates entropy of the input *x* using formula
%   XXX
%

%u = f./norm(f);
v = sort(abs(x(:)),'descend');

k = 1:length(v);
E = max((k.'.^(1/p)).*v);