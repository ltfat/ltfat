function [gamma,B,tt,cond0]=bspline_dual(a,alpha,beta,increase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       computes exponential B-spline B and dual                      %%%
%%%      window gamma for the Gabor frame G(B,alpha,beta)               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% a : vector of weights a_j = \lambda_j in Fourier representation of B
% alpha, beta : Gabor lattice parameters, must satisfy one of the cases:  
%               1.) 0 < alpha < length(a) and 0 < beta <= 1/length(a)
%               2.) alpha \in {1, ... , length(a)-1}, beta > 0, 
%                   alpha*beta < 1
%               3.) beta \in {1, 1/2, ... , 1/length(a)}, alpha > 0, 
%                   alpha*beta < 1
% increase : number of additional columns to compute window function
%            gamma; 'increased'-many are added at each side;
%            should be smaller than 100 to have comfortable execution-time
%
% OUTPUT:
% B : function values of EB-spline B
% gamma : each single row gamma(j,*) contains function values of a dual
%         Gabor window gamma, computed from pseudo-inverses of rectangular
%         blocks P0(x) of the matrix A(x);
%         P0(x) chooses columns kcol(j,1):kcol(j,2) of A(x), rows i1:i2 are
%         computed accordingly;
%         in everey row of gamma there were used 'increse' more columns of
%         A(x) to compute the new values of a dual-window according to B
% tt : abscissa for B and gamma (maybe only for g, if alpha*beta -> 0)
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
% discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = length(a);
x = ceil(max(1+(ceil(m/(1/beta-alpha))+increase+1)/alpha/beta-m/alpha,...
        (ceil(alpha/(1/beta-alpha))+increase+1)/alpha/beta))+1;
K = 100; % integer, for stepsize alpha/K of evaluation of B, gamma
t = 0:alpha/K:alpha*(1-1/K); % for stepping through the interval [0,alpha)
tt = -x*alpha:alpha/K:x*alpha; % choose same stepsize for t and tt
    % left and right bounds large enough for the support of gamma
gamma = zeros(1,length(tt)); % dual window
B = expbspval(a,tt); % EB-spline


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wait = waitbar(0,'Computing dual windows. Please wait...');

for s=1:length(t)
    % step through the interval [0,alpha)
    x0 = t(s); % compute dual window at points (x+j*alpha)
    
    % column indices for rectangular P0
    k = 1; l = 1;
    while x0+k*(1/beta-alpha) < m
        k = k + 1;
    end
    while x0-l*(1/beta-alpha) > 0
        l = l + 1;
    end
    k1 = k + increase;
    k2 = l + increase;
    
    % row indices for rectangular P0
    i1 = k1;
    i2 = k2;
    while x0-(i1+1)*alpha+(k1+1)/beta >= m
        i1 = i1 + 1;
    end
    while x0+(i2+1)*alpha-(k2+1)/beta <= 0
        i2 = i2 + 1;
    end
    
    % Computation of P0(x0)
    P0 = PreB(a,x0,alpha,beta,-i1,i2,-k1,k2);

    % computation of pseudo-inverse matrix of P0
    [P0inv,condP0] = my_pinv(P0);
    gamma(s-1+find(tt==0)-K*i1:K:s-1+find(tt==0)+K*i2) = beta*P0inv(k1+1,:); % row index
    % points to the "j=0" row of P0inv

    % record the condition number
    cond(s) = condP0;

    % actualizing waitbar
    counter = s;
    if ishandle(wait)
        waitbar(counter/K,wait,[num2str((counter*100)/K) '%']);
    else
        warndlg('Process was canceled by user.','!! Warning !!');
        return;
    end
end

cond0 = mean(cond);
% display('Minimal condition number of the matrices:')
% display(min(cond))
% display('Maximal condition number of the matrices:')
% display(max(cond))

try
    close(wait);
end

display('Test of Wexler-Raz-Biorthogonal condition:')
display((gamma*B')/beta/K)