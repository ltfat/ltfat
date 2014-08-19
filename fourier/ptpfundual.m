function [gamma,nlen]=ptpfundual(g,w,a,M,L,increase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       computes totally positive function g and dual                 %%%
%%%      window gamma for the Gabor frame G(g,alpha,beta)               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% w : vector of reciprocals a_j=1/delta_j in Fourier representation of g
% a : time shift, given by an integer number of sampling points
% M : number of channels
% L : length of a period
% increase : number of additional columns to compute window function
%            gamma; 'increased'-many are added at each side;
%            should be smaller than 100 to have comfortable execution-time
%
% OUTPUT:
% gamma : periodized dual window for the discrete TP function g with given
%         weihts w
%
%
%   References: grst13 kl12 bagrst14 klst14

% AUTHORS: Joachim Stoeckler, Tobias Kloos, 2012-2014




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chek input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = w/sqrt(L);
alpha = a;
beta = 1/M;
% check alpha beta
if (alpha<=0) || (beta<=0)
    error('lattice parameters alpha, beta must be positive')
end
if (alpha*beta>=1)
    error('lattice parameters must satisfy alpha*beta<1')
end

% compute m n and check that a has nonzero entries
w = sort(w(:)); % sort and make it a column vector
mult=myknt2mlt(w);
m = length(find(w>0));
n = length(find(w<0));
if (m+n)<length(w)
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
varl = floor((k1+m-1)/(alpha*beta))-1;
varr = ceil((k2-n+1)/(alpha*beta))+1;
x = varl*alpha:alpha:varr*alpha;
i0 = abs(varl)+1; % index of "central" row of P(x)
y = (k1-1)/beta:(1/beta):(k2+1)/beta;
k0 = abs(k1-1)+1; % index of "central" column of P(x)
[yy,xx] = meshgrid(y,x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:(a-1); % for stepping through the interval [0,alpha)
tt = varl*alpha:varr*alpha; % choose same stepsize for t and tt
    % left and right bounds large enough for the support of gamma
tt0 = abs(varl*a)+1; % index for tt == 0
gamma = zeros(1,length(tt)); % dual window


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k1 = k1+k0;
k2 = k2+k0;

for k=1:length(t)
    % step through the interval [0,alpha)
    x0 = t(k); % compute dual window at points (x+j*alpha)

    % row indices for rectangular P0
    if case0 == 0
      %  i1 = min(find((x+x0)>(y(k1+m-1))));
      i1 = floor((k1-k0+m-1)/alpha/beta-x0/alpha)+1+i0;
      %  i2 = max(find((x+x0)<(y(k2-n+1))));
      i2 = ceil((k2-k0-n+1)/alpha/beta-x0/alpha)-1+i0;
    else
      %  i1 = max(find((x+x0)<(y(k1))));
      i1 = ceil((k1-k0)/alpha/beta-x0/alpha)-1+i0;
      %  i2 = max(find((x+x0)<(y(k2))));
      i2 = ceil((k2-k0)/alpha/beta-x0/alpha)-1+i0;
    end

    % Computation of P0(x0)
    % z0 is the matrix of the abscissa x0+j*alpha-k/beta, j=i1:i2, k=k1:k2,
    % z1 puts all these abscissae into a row vector.
    % The computation of g(z1) is done as described above for the
    % vector tt.
    z0 = x0+xx(:,k1:k2)-yy(:,k1:k2);
    z1 = z0(:)';
    za = w*z1;     % matrix of values w(j)*z1(k)
    za = max(-1,za);  % for numerical stability / Inf could be created next step
    zb = exp(-za).*(za>=0).*repmat(sign(z1),length(w),1);
    % matrix of values for the divided difference
    m0 = find(z1==0); % compute values at 0 separately, one-sided only
    if ~isempty(m0)
        zb(:,m0)=repmat((w>0),1,length(m0));
    end
    if ~isempty(find(mult~=0)) % take the values of the derivates, if there are multiplicities
        col = find(mult~=0);
        for i = 1:length(col)
            zb(col(i),:)=((-z1).^(mult(col(i))).*exp(-za(col(i),:)).*(za(col(i),:)>=0).*sign(z1))./(factorial(mult(col(i))));
        end
    end
    c = (-1)^(m+n-1)*prod(w)*L^(1/4); % normalization constant for g
    A0 = c*divdiff_vector(w,zb);  % formula for g with divided differences
    A0 = reshape(A0,size(z0));
    P0 = A0(i1:i2,:);

    % computation of pseudo-inverse matrix of P0
    P0inv = pinv(P0);
    gamma(k-1+tt0-a*(i0-i1):a:k-1+tt0+a*(i2-i0)) = beta*P0inv(k0-k1+1,:); % row index k0-k1a+1
    % points to the "j=0" row of P0inv
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% periodization of gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nlen = length(gamma);
nr = ceil(nLen/L);
v = zeros(1,nr*L);
v(1:length(gamma)) = gamma;
v = [v(tt0:end),v(1:tt0-1)];

gamma = sum(reshape(v,L,nr),2);
gamma = gamma(:);

[scal,err]=gabdualnorm(g,gamma,a,M,L);
assert(err<1e-10,sprintf(['%s: Assertion failed. This is not a valid ',...
                          ' dual window.'],upper(mfilename)));

gamma = gamma/scal;

nlen = min([L,nlen]);

% TO DO: as gamma might have a compact support, nlen<L will indicate that.
% Howewer, since it is symmetric only in special cases, we cannot
% use middlepad(gamma,nlen) as in pbspline
%
% 2 options:
%
%    Either work with the offset as in struct('h',gamma,'offset',offset)
%     (this would require refactoring all gabor routines)
%
%    or pad the shorter "tail" to length of the longer tail.

