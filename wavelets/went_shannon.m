function E = went_shannon(x)
%WENT_SHANNON Shannon additive entropy
%   Usage:  E = went_shannon(x)
%
%   `E = went_shannon(x)` calculates entropy of the input *x* using formula
%   XXX
%
%u = x./norm(x);
u = x(x~=0).^2;
E = -sum(u.*log(u));