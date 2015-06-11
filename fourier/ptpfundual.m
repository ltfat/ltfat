function [gd,nlen]=ptpfundual(L,w,a,M,varargin)
%PTPFUNDUAL Sampled, periodized dual TP function of finite type
%   Usage: gd=ptpfundual(L,w,a,M)
%          gd=ptpfundual(L,w,a,M,inc)
%          gd=ptpfundual(L,w,a,M,...)
%          [gd,nlen]=ptpfundual(L,w,a,M,...)
%
%   Input parameters:
%         L     : Window length.
%         w     : Vector of reciprocals $w_j=1/\delta_j$ in Fourier representation of *g*
%         a     : Length of time shift.
%         M     : Number of channels.
%
%   Output parameters:
%         gd    : The periodized totally positive function dual window.
%         nlen  : Number of non-zero elements in gd.
%
%   `ptpfundual(L,w,a,M)` computes samples of a dual window of totally 
%   positive function of finite type >=2 with weights *w*. Please see
%   |ptpfun| for definition of totally positive functions.
%   The lattice parameters $a,M$ must satisfy $M \geq a$ to ensure the
%   system is a frame. 
%
%   `ptpfundual(L,w,a,M,inc)` works as above, but integer *inc* denotes
%   number of additional columns to compute window function *gd*. 
%   'inc'-many are added at each side. It should be smaller than 100 to
%   have comfortable execution-time. The higher the number the closer *gd*
%   is to the canonical dual window.
%   The default value is 0.
%
%   `[gd,nlen]=ptpfundual(...)` as *gd* migth have a compact support, 
%   *nlen* contain a number of non-zero elements in *gd*. This is the case 
%   when *gd* is symmetric. If *gd* is not symmetric, *nlen* is extended
%   to the twice the length of the the longer tail.   
%
%   If $nlen = L$, *gd* has a 'full' support meaning it is a periodization
%   of a dual TP function.
%
%   If $nlen < L$, additional zeros can be removed by calling 
%   `gd=middlepad(gd,nlen)`.
%   
%   In addition, `ptpfundual` accepts any of the flags from |normalize|. 
%   The output will be normalized as specified. The default normalization 
%   is `'energy'`.
%
%   If `ptpfundual` is intended to be used in conjuction with |ptpfun|,
%   an additional flag `matchscale` must be used. In such case, the 
%   original window is calculated and normalized according to the 
%   |normalize| flag and *gd* is scaled using |gabdualnorm|.
%
%   See also: dgt, ptpfun, gabdualnorm, normalize
%
%   References: grst13 kl12 bagrst14 klst14

% AUTHORS: Joachim Stoeckler, Tobias Kloos, 2012-2014

complainif_notenoughargs(nargin,4,upper(mfilename));
complainif_notposint(L,'L',upper(mfilename));
complainif_notposint(a,'a',upper(mfilename));
complainif_notposint(M,'M',upper(mfilename));

% Check lattice
if M<=a
    error('%s: Lattice parameters must satisfy M>a.',upper(mfilename));
end

% Check tpfun vect.
if isempty(w) || ~isnumeric(w)
    error('%s: w must be a nonempty numeric vector.', upper(mfilename));
end

if numel(w)<2
    error(['%s: The tp fun. finite type must be >=2 (number of ',...
           'elements of w).'], upper(mfilename));
end

if any(w==0)
    error('%s: All weights w must be nonzero.', upper(mfilename));
    % TO DO: Also add a warning if w is very small or big?
end

% Define initial value for flags and key/value pairs.
definput.import={'normalize'};
definput.importdefaults={'null'};
definput.keyvals.inc = 0;
definput.flags.scale = {'nomatchscale','matchscale'};
[flags,~,inc]=ltfatarghelper({'inc'},definput,varargin);

complainif_notnonnegint(inc,'inc',upper(mfilename));

% TP functions are scale invariant so we do scaling directly on w.
wloc = w/sqrt(L);
% Converting a, M to alpha, beta
alpha = a;
beta = 1/M;

% compute m n and check that a has nonzero entries
wloc = sort(wloc(:)); % sort and make it a column vector
mult=wmult(wloc);
m = length(find(wloc>0));
n = length(find(wloc<0));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations specially for computation of gd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = floor(1/(1-alpha*beta)+eps);

% check special cases according to n and m
case0 = 0;
if m == 0
    if n >= 2
        N = (n-1)*(r+1); % minimal column size
        k1 = 0;
        k2 = N-1;
        case0 = 0;
    end
elseif n == 0
    if m >= 2
        N = (m-1)*(r+1); % minimal column size
        k2 = 0;
        k1 = -(N-1);
        case0 = 0;
    end
else
    k1 = -m*(r+1)+1;     % column index k1 from the paper
    k2 = n*(r+1)-1;      % symmetric choice
    N = k2-k1+1;
end

k1 = k1-inc;
k2 = k2+inc;

% minimal values for x and y
varl = floor((k1+m-1)/(alpha*beta))-1;
varr = ceil((k2-n+1)/(alpha*beta))+1;
x = linspace(varl*alpha,varr*alpha,varr-varl+1);
i0 = abs(varl)+1; % index of "central" row of P(x)
y = linspace((k1-1)*M,(k2+1)*M,k2-k1+3);
k0 = abs(k1-1)+1; % index of "central" column of P(x)
[yy,xx] = meshgrid(y,x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:(a-1); % for stepping through the interval [0,alpha)
tt = varl*alpha:varr*alpha; % choose same stepsize for t and tt
    % left and right bounds large enough for the support of gamma
tt0 = abs(varl*a)+1; % index for tt == 0
gd = zeros(1,length(tt)); % dual window


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
    elseif case0 == 1  % m=0
      %  i1 = max(find((x+x0)<(y(k1))));
      i1 = ceil((k1-k0)/alpha/beta-x0/alpha)-1+i0;
      %  i2 = max(find((x+x0)<(y(k2-n+1))));
      i2 = ceil((k2-k0-n+1)/alpha/beta-x0/alpha)-1+i0;
    elseif case0 == 2  % n=0
      %  i1 = min(find((x+x0)>(y(k1+m-1))));
      i1 = floor((k1-k0+m-1)/alpha/beta-x0/alpha)+1+i0;
      %  i2 = min(find((x+x0)>(y(k2))));
      i2 = floor((k2-k0)/alpha/beta-x0/alpha)+1+i0;
    end

    % Computation of P0(x0)
    % z0 is the matrix of the abscissa x0+j*alpha-k/beta, j=i1:i2, k=k1:k2,
    % z1 puts all these abscissae into a row vector.
    % The computation of g(z1) is done as described above for the
    % vector tt.
    z0 = x0+xx(:,k1:k2)-yy(:,k1:k2);
    z1 = z0(:)';
    za = wloc*z1;     % matrix of values w(j)*z1(k)
    za = max(-1,za);  % for numerical stability / Inf could be created next step
    zb = exp(-za).*(za>=0).*repmat(sign(z1),length(wloc),1);
    % matrix of values for the divided difference
    m0 = find(z1==0); % compute values at 0 separately, one-sided only
    if ~isempty(m0)
        if m>0
            zb(:,m0)=repmat((wloc>0),1,length(m0));
        else
            zb(:,m0)=repmat((wloc<0),1,length(m0));
        end
    end
    if ~isempty(find(mult~=0)) % take the values of the derivates, if there are multiplicities
        col = find(mult~=0);
        for i = 1:length(col)
            zb(col(i),:)=((-z1).^(mult(col(i))).*zb(col(i),:))./(factorial(mult(col(i))));
        end
    end
    c = (-1)^(m+n-1)*prod(wloc)*L^(1/4); % normalization constant for g
    A0 = divdiff_vector(wloc,mult,zb);  % formula for g with divided differences
    A0 = reshape(A0,size(z0));
    P0 = A0(i1:i2,:);
    % computation of pseudo-inverse matrix of P0
    P0inv = pinv(P0);
    gd(k-1+tt0-a*(i0-i1):a:k-1+tt0+a*(i2-i0)) = beta*P0inv(k0-k1+1,:); % row index k0-k1a+1
    % points to the "j=0" row of P0inv
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% periodization of gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nlen = length(gd);
nr = ceil(nlen/L);
v = zeros(1,nr*L);
v(1:length(gd)) = gd;
v = [v(tt0:end),v(1:tt0-1)];

gd = sum(reshape(v,L,nr),2);
gd = gd(:);

% Determine nlen
if nlen<L
   negsupp = tt0-1;
   possupp = nlen-tt0; % excluding zero pos.
   nlen = 2*max([negsupp,possupp])+1;
end
nlen = min([L,nlen]);

if flags.do_matchscale
   g = ptpfun2(L,w,flags.norm);
   [scal,err]=gabdualnorm(g,gd,a,M,L);
   assert(err<1e-10,sprintf(['%s: Assertion failed. This is not a valid ',...
                             ' dual window.'],upper(mfilename)));
   gd = gd/scal;
else
   gd = normalize(gd,flags.norm); 
end


