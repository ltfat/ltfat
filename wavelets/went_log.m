function E = went_log(x)
%WENT_LOG Logarithmic aditive entropy
%   Usage:  E = went_log(x)
%
%   `E = went_log(x)` calculates entropy of the input *x* using formula
%   XXX
%

%u = x(:)./norm(x(:));
%u = u(u~=0);
x = x(x~=0);
E = sum(log(x(:).^2));