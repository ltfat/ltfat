function [g,gamma,tt,cond0,wnorm]=TPframe(a,alpha,beta,increase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       computes totally positive function g and dual                 %%%
%%%      window gamma for the Gabor frame G(g,alpha,beta)               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% a : vector of reciprocals a_j=1/delta_j in Fourier representation of g
% alpha, beta : Gabor lattice parameters, must satisfy 0 < alpha,beta and
%               alpha*beta < 1
% increase : number of additional columns to compute window function
%            gamma; 'increased'-many are added at each side;
%            should be smaller than 100 to have comfortable execution-time
%
% OUTPUT:
% g : function values of TP function g
% gamma : each single row gamma(j,*) contains function values of a dual
%         Gabor window gamma, computed from pseudo-inverses of rectangular
%         blocks P0(x) of the matrix A(x);
%         P0(x) chooses columns kcol(j,1):kcol(j,2) of A(x), rows i1:i2 are
%         computed accordingly;
%         in everey row of gamma there were used 'increse' more columns of
%         A(x) to compute the new values of a dual-window according to g
% tt : abscissa for g and gamma (maybe only for g, if alpha*beta -> 0)
% cond0: gives max. condition number of P0(x)
% wnorm : wiener norm of gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   References:
%     K. Groechenig and J. Stoeckler,
%     Gabor Frames and Totally Positive Functions,
%     Duke Mathematical Journal, Volume 162, pp.1003-1031, 2013.
%
%     T. Kloos,
%     Gabor Frames total-positiver Funktionen endlicher Ordnung,
%     Diploma thesis, University of Dortmund, 2012.
%
%     S. Bannert and K. Groechenig and J. Stoeckler,
%     Discretized Gabor Frames of Totally Positive Functions,
%     IEEE Transactions on Information Theory, Volume 60, pp.159-169, 2014.
%
%     T. Kloos and J. Stoeckler,
%     Zak transforms and Gabor frames of totally positive functions and 
%     exponential B-splines,
%     Journal of Approximation Theory, Volume 184, pp.209-237, 2014.

%   (c) Joachim Stoeckler,
%       Tobias Kloos, 2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chek input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check alpha beta
if (alpha<=0) || (beta<=0)
    error('lattice parameters alpha, beta must be positive')
end
if (alpha*beta>=1)
    error('lattice parameters must satisfy alpha*beta<1')
end

% compute m n and check that a has nonzero entries
a = sort(a(:)); % sort and make it a column vector
mult=myknt2mlt(a);
m = length(find(a>0));
n = length(find(a<0));
if (m+n)<length(a)
    error('all entries of the vector a must be non-zero')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations specially for computation of gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = floor(1/(1-alpha*beta)+eps);

% check special cases according to n and m
case0 = 0;
if m == 0
    if n >= 2
        N = (n-1)*(r+1); % minimal column size
        k1 = 0;
        k2 = N-1;
        case0 = 1;
    else
        error('g should be a totally positive function of finite type >= 2')
    end
elseif n == 0
    if m >= 2
        N = (m-1)*(r+1); % minimal column size
        k1 = 0;
        k2 = N-1;
        case0 = 1;
    else
        error('g should be a totally positive function of finite type >= 2')
    end
elseif n == 1
    N = m*(r+1)+1;       % minimal column size
    k1 = -m*(r+1)+1;     % column index k1 from the paper
    k2 = k1+N-1;         % column index k2 from the paper
else
    N = (m+n-1)*(r+1);   % minimal column size
    k1 = -m*(r+1)+1;     % column index k1 from the paper
    k2 = k1+N-1;         % column index k2 from the paper
end

k1 = k1-increase;
k2 = k2+increase;

% minimal values for x and y
var = max(abs(floor((k1+m-1)/(alpha*beta))),abs(ceil((k2-n+1)/(alpha*beta))));
x = -var*alpha:alpha:var*alpha;
i0 = min(find(x>=0)); % index of "central" row of P(x)
y = (k1-1)/beta:(1/beta):(k2+1)/beta;
k0 = min(find(y>=0)); % index of "central" column of P(x)
[yy,xx] = meshgrid(y,x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = 100; % integer, for stepsize alpha/K of evaluation of g, gamma
t = 0:alpha/K:alpha*(1-1/K); % for stepping through the interval [0,alpha)
tt = min(-var,-100)*alpha:alpha/K:max(var,100)*alpha; % choose same stepsize for t and tt
    % left and right bounds large enough for the support of gamma
gamma = zeros(1,length(tt)); % dual window


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of g at the abscissae tt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The values g(tt(j)) are computed by the explicit formula in my note using 
% divided differences.
% Starting values for the divided differences are stored in matrix zb,
% the divided differences are computed with the usual recursive formula
% acting on each column zb(:,j) in order to compute g(tt(j)).
% This works only for pairwise distinct values in vector a, 
% it should be extended to higher multiplicity later.
za = a*tt;     % matrix of values a(j)*tt(k)
zb = exp(-za).*(za>=0).*repmat(sign(tt),length(a),1);
    % matrix of values for the divided difference
m0 = find(tt==0); % compute values at 0 separately, one-sided only
if ~isempty(m0)
    zb(:,m0)=repmat((a>0),1,length(m0));
end
if ~isempty(find(mult~=0)) % take the values of the derivates, if there are multiplicities
   col = find(mult~=0);
   for i = 1:length(col)
       zb(col(i),:)=((-tt).^(mult(col(i))).*exp(-za(col(i),:)).*(za(col(i),:)>=0).*sign(tt))./(factorial(mult(col(i))));
   end
end
c = (-1)^(m+n-1)*prod(a); % normalization constant for g
g = c*divdiff_vector(a,zb);  % the new formula for g in my note is used here


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wait = waitbar(0,'Computing dual windows. Please wait...');
k1 = k1+k0;
k2 = k2+k0;

for k=1:length(t)
    % step through the interval [0,alpha)
    x0 = t(k); % compute dual window at points (x+j*alpha)

    % row indices for rectangular P0
    if case0 == 0
        i1 = min(find((x+x0)>(y(k1+m-1))));
        i2 = max(find((x+x0)<(y(k2-n+1))));
    else
        i1 = max(find((x+x0)<(y(k1))));
        i2 = max(find((x+x0)<(y(k2))));
    end

    % Computation of P0(x0)
    % z0 is the matrix of the abscissa x0+j*alpha-k/beta, j=i1:i2, k=k1:k2,
    % z1 puts all these abscissae into a row vector.
    % The computation of g(z1) is done as described above for the
    % vector tt.
    z0 = x0+xx(:,k1:k2)-yy(:,k1:k2);
    z1 = z0(:)';
    za = a*z1;     % matrix of values a(j)*z1(k)
    za = max(-1,za);  % for numerical stability / Inf could be created next step
    zb = exp(-za).*(za>=0).*repmat(sign(z1),length(a),1);
    % matrix of values for the divided difference
    m0 = find(z1==0); % compute values at 0 separately, one-sided only
    if ~isempty(m0)
        zb(:,m0)=repmat((a>0),1,length(m0));
    end
    if ~isempty(find(mult~=0)) % take the values of the derivates, if there are multiplicities
        col = find(mult~=0);
        for i = 1:length(col)
            zb(col(i),:)=((-z1).^(mult(col(i))).*exp(-za(col(i),:)).*(za(col(i),:)>=0).*sign(z1))./(factorial(mult(col(i))));
        end
    end
    A0 = c*divdiff_vector(a,zb);  % formula for g with divided differences
    A0 = reshape(A0,size(z0));
    P0 = A0(i1:i2,:);

    % computation of pseudo-inverse matrix of P0
    [P0inv,condP0] = my_pinv(P0);
    gamma(k-1+find(tt==0)-K*(i0-i1):K:k-1+find(tt==0)+K*(i2-i0)) = beta*P0inv(k0-k1+1,:); % row index k0-k1a+1
    % points to the "j=0" row of P0inv
    % record the condition number
    cond(k) = condP0;

    % actualizing waitbar
    counter = k;
    if ishandle(wait)
        waitbar(counter/K,wait,[num2str((counter*100)/K) '%']);
    else
        warndlg('Process was canceled by user.','!! Warning !!');
        return;
    end
end

cond0 = mean(cond);

try
    close(wait);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of the wiener norm and framebounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wnorm = 0;
for i = 1:K:(length(gamma(:))-K)
    wnorm = wnorm + max(abs(gamma(i:i+K)));
end

display('Test of Wexler-Raz-Biorthogonal condition:')
display((gamma*g')/(beta*100))