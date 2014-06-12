function [P,P0] = PreB(a,x0,alpha,beta,r1,r2,c1,c2)

% computes the submatrix P0(x0) with rows r1-r2 and columns c1-c2 
% of pre-gramian matrix P(x0) of the 
% exponential b-spline with vector of knots a

%   References:
%     O. Christensen and P. Massopust,
%     Exponential B-splines and the partition of unity property,
%     Adv. Comput. Math., Volume 37, pp.301-318, 2012.

%   (c) Tobias Kloos, 2013

if(nargin < 5)
   m0 = ceil(x0/(1/beta-alpha))-1;
   r1 = 0; r2 = m0;
   c1 = -1; c2 = m0+1;
end

%m0 = ceil(x0/(1/beta-alpha))-1
%x0/(1/beta-alpha)

% m = length(a);
% s = ceil(1/beta/alpha-1);
% l0 = ceil(x0/alpha); xo = x0-(l0-1)*alpha;
% k0 = ceil((m-x0)/alpha); xu = x0+(k0-1)*alpha;
% l1 = ceil((m-xo)/(1/beta-s*alpha));
% k1 = ceil(xu/(1/beta-s*alpha));

k = [r1:r2]; j = [c1:c2];
% k = [-l1*s-l0+1:k1*s+k0-1]; j = [-l1+1:k1-1];

P0 = x0 + repmat(alpha*k',1,length(j)) - repmat(j./beta,length(k),1);
P = expbspval(a,P0);