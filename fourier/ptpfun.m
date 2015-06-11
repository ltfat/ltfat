function g = ptpfun(L,w,varargin)
%PTPFUN Sampled, periodized totally positive function of finite type
%   Usage: g=ptpfun(L,w)
%          g=ptpfun(L,w,...)
%
%   Input parameters:
%         L    : Window length.
%         w    : Vector of reciprocals $w_j=1/\delta_j$ in Fourier representation of *g*
%
%   Output parameters:
%         g    : The periodized totally positive function.
%
%   `ptpfun(L,w)` computes samples of a periodized totally positive
%   function of finite type >=2 with weights *w* such that the Fourier
%   representation of the continuous TP function is given as:
%
%   ..             m 
%      ghat(f) = prod (1+2pijf/w(i))^(-1), 
%                 i=1  
%
%   .. math:: \hat{g}(\xi)=\prod_{i=1}^{m}\left(1+2\pi i j\xi /w(i)\right)^{-1},
%
%   where $m$=`numel(w)`$\geq 2$. The samples are the values of
%   the Zak transform Zg(x,0) with frequency variable set to zero. 
%   
%   *w* controls the function decay in the time domain. More specifically
%   the function decays as $exp(-min(wpos)*x)$ for $x->\infty$ 
%   and $exp(max(wneg)*x)$ for $x->-\infty$ assuming *w* contains 
%   positive numbers (wpos) and negative numbers (wneg).
%   If *w* has only positive (negative) values, the continuous 
%   TP function is zero on the negative (positive) axis.
%
%   In addition `ptpfun` accepts any of the flags from |normalize|. The 
%   output will be normalized as specified. The default normalization flag
%   is `'null'` i.e. no normalization.
%
%   See also: dgt, ptpfundual, gabdualnorm, normalize
%
%   References: grst13 kl12 bagrst14 klst14
%

%   AUTHORS: Joachim Stoeckler, Tobias Kloos  2012, 2014

complainif_notenoughargs(nargin,2,upper(mfilename));
complainif_notposint(L,'L',upper(mfilename));

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
[flags,keyvals]=ltfatarghelper({},definput,varargin);

w = sort(w(:))*sqrt(L);
m = length(w);
[wm,w] = wmult(w);
wmax = max(wm)+1;   % maximal multiplicity in w

% if k is the multiplicity of y in w, then we need derivatives
% up to order k-1 of the function
%    f(y) = exp(-x*y)/(1-exp(-y)) 
% where x is a parameter in [0,1].
% We use the Leibniz formula and explicit representations of 
% the derivatives of f1(y)=(1-exp(-y))^(-1) in terms of Euler-Frobenius
% polynomials:
% (d^(k-1))/(dy^(k-1)) f1(y) = (-1)^(k-1)*(1-exp(-y))^(-k)*sum_(j=1)^(k) pw(k,j)*exp(-(k-j)*y)
pw = zeros(wmax);   % square matrix of Euler-Frobenius coefficients 
pw(1,1) = 1;
for k = 2:wmax
   pw(k,1:k)=[k-1:-1:0].*[0,pw(k-1,1:k-1)]+[1:k].*[pw(k-1,1:k-1),0];
end  
% stable computation for x in [0,1] and all real y=w(j) of 
%     exp(-x*y)*(1-exp(-y))^(-k)*sum_(j=1)^(k) pw(k,j)*exp(-(k-j)*y)
x=[0:L-1]/L;
f1=zeros(m,L);
h=zeros(wmax,L);
for r=1:m
    wr = w(r);
    k = wm(r);
    if k==0
        if wr>0
            h(1,:) = exp(-wr*x);
            f1(r,:) = 1/(1-exp(-wr))*h(1,:);
        else
            h(1,:) = exp(wr*(1-x));
            f1(r,:) = 1/(exp(wr)-1)*h(1,:);
        end
    else  % only for multiple wr
        h(k+1,:) = exp(-abs(wr))*h(k,:);   
        if wr>0
            f1(r,:) = ((1-exp(-wr))^(-k-1)*pw(k+1,k+1:-1:1))*h(1:k+1,:);
        else
            f1(r,:) = ((exp(wr)-1)^(-k-1)*pw(k+1,1:k+1))*h(1:k+1,:);
        end
    end
end
% powers of x 
f2=ones(wmax,L);
for r=2:wmax
    f2(r,:)=x.*f2(r-1,:);
end
% input for divided difference
y=zeros(m,L);
for r=1:m
    k = wm(r);
    for ell=0:k    % Leibniz rule
        y(r,:) = y(r,:)+nchoosek(k,ell)*f1(r-k+ell,:).*f2(k-ell+1,:);
    end
    y(r,:) = (-1)^(k)/factorial(k)*y(r,:);
end

g = (-1)^(m-1) * prod(w) * divdiff_vector(w,wm,y) * L^(-3/4);
g = g(:);
g = normalize(g(:),flags.norm);



