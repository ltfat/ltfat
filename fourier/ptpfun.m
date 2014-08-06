function g = ptpfun(L,w,varargin)
%PTPFUN Sampled, periodized totally positive function of finite type
%   Usage: g=ptpfun(L,w)
%          g=ptpfun(L,w,...)
%
%   Input parameters:
%         L    : Length of vector.
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
%   where $m$=`numel(w)`$\geq 2$. The samples are obtained from the Zak
%   transform of the function.
%   
%   *w* controls the function decay in the time domain. More specifically
%   the function decays as exp(max(w)x) for x->\infty and exp(min(w)x) for
%   x->-\infty assuming w contains both positive and negative numbers.
%
%   In addition `ptpfun` accepts any of the flags from |normalize|. The 
%   output will be normalized as specified.
%
%   See also: dgt
%
%   References: grst13 kl12 bagrst14 klst14
%

%   AUTHORS: Joachim Stoeckler, Tobias Kloos  2012, 2014

complainif_notenoughargs(nargin,2,upper(mfilename));

if numel(w)<2
    error(['%s: The tp fun. finite type must be >=2 (number of ',...
           'elements of w).'], upper(mfilename));
end

if any(w==0)
    error('%s: All weights w must be nonzero.', upper(mfilename));
end

if all(w<0) || all(w>0)
    error(['%s: Only positive or only negative weights w are not ',...
           'supported yet.'], upper(mfilename));
end


% Define initial value for flags and key/value pairs.
definput.import={'normalize'};
[flags,keyvals]=ltfatarghelper({},definput,varargin);

w = sort(w(:))*sqrt(L);
m = length(w);
y = w*[0:L-1]/L;
x = exp(-y) ./ repmat(1-exp(-w),1,L);
g = (-1)^(m-1) * prod(w) * divdiff_vector(w,x) * L^(-3/4);

g = normalize(g(:),flags.norm);


