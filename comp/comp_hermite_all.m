function y = comp_hermite_all(n, x)
%COMP_HERMITE_ALL  Compute all Hermite functions up to an order
%   Usage:  y = hermite_fun_all(n, x);
%
%   This function evaluates the Hermite functions
%   of degree 0 through n-1 at the vector x.
%   The functions are normalized to have the L^2 norm
%   on (-inf,inf) equal to one. No effort is made to 
%   avoid unerflow during recursion.	
%   
%   Input parameters: 
%     n     : the number of Hermite functions
%     x     : the vector of arguments
%
%   Output parameters:
%     y     : the values of the first n Hermite functions at 
%             the nodes x


%   T. Hrycak, Mar. 22, 2006
%   Last modified   July 17, 2007
%
%
%       

rt = 1 / sqrt(sqrt(pi));

%     
%     conducting the recursion.
%

y = zeros(length(x), n);

y(:, 1) = rt * exp(-0.5 * x.^2);
if n > 1
   y(:, 2) = rt * sqrt(2) * x .* exp(-0.5 * x.^2);
end
for k = 2:n-1
        y(:, k+1) = sqrt(2)*x.*y(:, k) - sqrt(k-1)*y(:, k-1);
        y(:, k+1) = y(:, k+1)/sqrt(k);
end


