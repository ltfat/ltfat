function gamma = pebfundual(w,a,M,L,width,increase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PEBFUNDUAL Dual window of sampled, periodized EB-spline
%   Usage: g=pebfundual(w,a,M,L,width,increase)
%          g=pebfundual(w,a,M,L,width,...)
%          g=pebfundual(w,a,M,L,...)
%
%
% INPUT:
% w : vector of weights of g
% a : time shift, given by an integer number of sampling points
% M : number of channels
% L : length of a period
% width: integer stretching factor of the window *g*
%
% increase : number of additional columns to compute window function
%            gamma; 'increased'-many are added at each side;
%            should be smaller than 100 to have comfortable execution-time
%
% OUTPUT:
% gamma : periodized dual window for the discrete EB-spline g with given
%         weights w
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
%
%     T. Kloos and J. Stoeckler and K. Groechenig,
%     Implementation of discretized Gabor frames and their duals
%     IEEE Transactions on Information Theory, 2016.

%   (c) Joachim Stoeckler,
%       Tobias Kloos, 2012-2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(w) || ~isnumeric(w)
    error('%s: w must be a nonempty numeric vector.', upper(mfilename));
end

if nargin == 4
    width = floor(sqrt(L));
    increase = 10;
end
if nargin == 5
    increase = 10;
end

w = sort(w(:)); % sort and make it a column vector
m = length(w);
alpha = a/width;
beta = width/M;

% check alpha beta
if (alpha<=0) || (beta<=0)
    error('lattice parameters alpha, beta must be positive')
end
if (alpha*beta>=1)
    error('lattice parameters must satisfy alpha*beta<1')
end
if (alpha >= m)
    error('a/width must be smaller than length(w)')
end
check = 1;
if (1/beta == floor(1/beta))
    check = 0;
elseif (alpha == floor(alpha))
    check = 0;
elseif (1/alpha == floor(1/alpha)) && (beta < 1)
    check = 0;
end
if check == 1
    warning('output may not be a dual window; M/width should be a small integer')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations specially for computation of gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1 = -(ceil(beta*(1/2*(m+alpha)+alpha*increase)/(1-alpha*beta))-1);
k2 = ceil(beta*(1/2*(m+alpha)+alpha*increase)/(1-alpha*beta))-1;
i1 = floor(m/2/alpha+(k1-1)/(alpha*beta)-1);
i2 = ceil((k2+1)/(alpha*beta)-(m-alpha)/2/alpha+1);

% minimal values for x and y
x = (i1-1)*alpha:alpha:(i2+1)*alpha;
i0 = abs(i1-1)+1; % index of "central" row of P(x)
y = k1/beta:1/beta:k2/beta;
k0 = abs(k1)+1; % index of "central" column of P(x)
[yy,xx] = meshgrid(y,x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:1/width:(a-1)/width;% for stepping through the interval [0,alpha)
j0 = ceil((m*width-a)/2);
t = t + j0/width;
tt = i1*alpha:1/width:i2*alpha; % choose same stepsize for t and tt
% left and right bounds large enough for the support of gamma
tt0 = abs(i1*a)+1; % index for tt == 0
gamma = zeros(1,length(tt)); % dual window
c0 = zeros(1,length(t));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%wait = waitbar(0,'Computing dual windows. Please wait...');

for k = 1:length(t)
    % step through the interval [0,alpha)
    x0 = t(k); % compute dual window at points (x+j*alpha)
    
    k1x = -(ceil((beta*(m-x0+alpha*increase))/(1-alpha*beta))-1) + k0;
    k2x = ceil((beta*(x0+alpha*increase))/(1-alpha*beta))-1 + k0;
    
    i1x = ceil((m*beta+(k1x-k0)-1-x0*beta)/(alpha*beta)) + i0;
    i2x = floor(((k2x-k0)+1-x0*beta)/(alpha*beta)) + i0;
    
    
    % Computation of P0(x0)
    % z0 is the matrix of the abscissa x0+j*alpha-k/beta, j=i1:i2, k=k1:k2,
    % z1 puts all these abscissae into a row vector.
    % The computation of g(z1) is done as described above for the
    % vector tt.
    z0 = x0+xx(:,k1x:k2x)-yy(:,k1x:k2x);
    z1 = z0(:)';
    
    lz1 = length(z1);
    z1 = repmat(z1.',1,m) + repmat([-m+1:0],lz1,1);
    z1 = z1(:)';
    Y = zeros(m-1,length(z1));
    
    for q = 1:m-1
        if w(q) == w(q+1)
            Y(q,:) = z1.*exp(w(q)*z1).*(z1>=0).*(z1<=1) + ...
                (2-z1).*exp(w(q)*z1).*(z1>1).*(z1<=2);
        else
            Y(q,:) = (exp(w(q)*z1)-exp(w(q+1)*z1))/(w(q)-w(q+1)).*(z1>=0).*(z1<=1) + ...
                (exp(w(q)-w(q+1))*exp(w(q+1)*z1)-exp(w(q+1)-w(q))*exp(w(q)*z1))/(w(q)-w(q+1)).*(z1>1).*(z1<=2);
        end
    end
    
    for q = 2:m-1
        for j = 1:m-q
            if w(j) == w(j+q)
                Y(j,(q-1)*lz1+1:end) = z1((q-1)*lz1+1:end)/q .* Y(j,(q-1)*lz1+1:end) + ...
                    exp(w(j))*(q+1-z1((q-1)*lz1+1:end))/q .* Y(j,(q-2)*lz1+1:end-lz1);
            else
                Y(j,(q-1)*lz1+1:end) = ( Y(j,(q-1)*lz1+1:end) - Y(j+1,(q-1)*lz1+1:end) + ...
                    exp(w(j))*Y(j+1,(q-2)*lz1+1:end-lz1) - exp(w(j+q))*Y(j,(q-2)*lz1+1:end-lz1) )/ ...
                    (w(j)-w(j+q));
            end
        end
    end
    
    if m == 1
        A0 = exp(w*z1).*(z1>=0).*(z1<=1);
    else
        A0 = Y(1,end-lz1+1:end);
    end
    
    A0 = reshape(A0,size(z0));
    P0 = A0(i1x:i2x,:);
    
    % computation of pseudo-inverse matrix of P0
    P0inv = pinv(P0);
    
    gamma(k-1+tt0+j0+a*(i1x-i0):a:k-1+tt0+j0+a*(i2x-i0)) = beta*P0inv(k0-k1x+1,:); % row index k0-k1a+1
    % points to the "j=0" row of P0inv
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% periodization of gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nr = ceil(length(gamma)/L);
v = zeros(1,nr*L);
v(1:length(gamma)) = gamma;
v = [v(tt0:end),v(1:tt0-1)];

gamma = sum(reshape(v,L,nr)',1);
gamma = gamma/sqrt(width);