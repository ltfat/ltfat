function g = pebfun(w,L,width)
%PEBFUN Sampled, periodized EB-spline
%   Usage: g=pebfun(w,L,width)
%          g=pebfun(w,L,...)
%
%   Input parameters:
%         L    : Window length.
%         w    : Vector of weights of *g*
%         width: integer stretching factor, the support of g is width*length(w)
%
%   Output parameters:
%         g    : The periodized EB-spline.
%
%   See also: dgt, pebfundual, gabdualnorm, normalize
%
%   References: grst13 kl12 bagrst14 klst14
%

%   AUTHORS: Joachim Stoeckler, Tobias Kloos  2012, 2016

if isempty(w) || ~isnumeric(w)
    error('%s: w must be a nonempty numeric vector.', upper(mfilename));
end

w = sort(w);
if nargin == 2
    width = floor(sqrt(L));
end
n = length(w);
x = linspace(0,n-1/width,n*width);
x = x(:);
m = length(x);
x = repmat(x,1,n) + repmat([-n+1:0],m,1);
x = x(:)';
Y = zeros(n-1,length(x));

for k = 1:n-1
    if w(k) == w(k+1)
        Y(k,:) = x.*exp(w(k)*x).*(x>=0).*(x<=1) + ...
            (2-x).*exp(w(k)*x).*(x>1).*(x<=2);
    else
        Y(k,:) = (exp(w(k)*x)-exp(w(k+1)*x))/(w(k)-w(k+1)).*(x>=0).*(x<=1) + ...
            (exp(w(k)-w(k+1))*exp(w(k+1)*x)-exp(w(k+1)-w(k))*exp(w(k)*x))/(w(k)-w(k+1)).*(x>1).*(x<=2);
    end
end

for k = 2:n-1
    for j = 1:n-k
        if w(j) == w(j+k)
            Y(j,(k-1)*m+1:end) = x((k-1)*m+1:end)/k .* Y(j,(k-1)*m+1:end) + ...
                exp(w(j))*(k+1-x((k-1)*m+1:end))/k .* Y(j,(k-2)*m+1:end-m);
        else
            Y(j,(k-1)*m+1:end) = ( Y(j,(k-1)*m+1:end) - Y(j+1,(k-1)*m+1:end) + ...
                exp(w(j))*Y(j+1,(k-2)*m+1:end-m) - exp(w(j+k))*Y(j,(k-2)*m+1:end-m) )/ ...
                (w(j)-w(j+k));
        end
    end
end

if n == 1
    y = exp(w*x).*(x>=0).*(x<=1);
else
    y = Y(1,end-m+1:end);
end

if m <= L
    g = [y,zeros(1,L-m)]/sqrt(width);
else
    y = [y,zeros(1,ceil(m/L)*L-m)];
    g = sum(reshape(y,L,ceil(m/L)),2).'/sqrt(width);
end

end