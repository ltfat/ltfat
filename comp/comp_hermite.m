function y = comp_hermite(n, x);
%COMP_HERMITE   Compute sampling of continuous Hermite function.
%   Usage:  y = comp_hermite(n, x);
%
%   COMP_HERMITE(n, x) evaluates the n-th Hermite function at the vector x.
%   The function is normalized to have the L^2(-inf,inf) norm equal to one.
%
%   A minimal effort is made to avoid underflow in recursion.
%   If used to evaluate the Hermite quadratures, it works for n <= 2400
%

% AUTHOR:
%   T. Hrycak, Oct 5, 2005
%   Last modified July 17, 2007


rt = 1 / sqrt(sqrt(pi));

if n == 0
  y = rt * exp(-0.5 * x.^2);
end
if n == 1
  y = rt * sqrt(2) * x .* exp(-0.5 * x.^2);
end

%     
%     if n > 2, conducting the recursion.
%

if n >= 2
  ef = exp(-0.5 * (x.^2) / (n+1));
  tmp1 = rt * ef;
  tmp2 = rt * sqrt(2) * x .* (ef.^2);
  for k = 2:n
    y = sqrt(2)*x.*tmp2 - sqrt(k-1)*tmp1 .* ef;
    y = ef .* y / sqrt(k);
    tmp1 = tmp2;
    tmp2 = y;
  end
end
  




