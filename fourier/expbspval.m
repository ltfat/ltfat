function y = expbspval(a,x)
% computation of EB-splines with distinct weights

%   References:
%     O. Christensen and P. Massopust,
%     Exponential B-splines and the partition of unity property,
%     Adv. Comput. Math., Volume 37, pp.301-318, 2012.

%   (c) Tobias Kloos, 2013

if sum(knt2mlt(a)) > 0
    error('all knots must be distinct')
end

a = sort(a);
m = length(a);
al = 1:m;
y = zeros(size(x));

for k = 1:m
    xk = find(x >= k-1 & x < k);
    if isempty(xk) == 0;
        comb = combine(m,k-1);
        for l = 1:m
            if k ~= 1
                s = sum(prod(reshape(exp(a(comb(sum(comb==l,2)==0,:))),size(comb(sum(comb==l,2)==0,:))),2));
            else
                s = 1;
            end
            p = prod((a(l)-a(find(al ~= l))));
            y(xk) = y(xk) + (s/p * exp(a(l)*(x(xk)-k+1)));
        end
    end
    y(xk) = (-1)^(k-1) * y(xk);
end

end