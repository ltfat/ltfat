function c=filterbank(f,g,a,varargin);  
%FILTERBANK   Apply filterbank
%   Usage:  c=filterbank(f,g,a);
%
%   `filterbank(f,g,a)` applies the filters given in *g* to the signal
%   *f*. Each subband will be subsampled by a factor of *a* (the
%   hop-size). In contrast to |ufilterbank|_, *a* can be a vector so the
%   hop-size can be channel-dependant. If *f* is a matrix, the
%   transformation is applied to each column.
%
%   The filters *g* must be a cell-array, where each entry in the cell
%   array corresponds to an FIR filter.
%
%   The output coefficients are stored a cell array. More precisely, the
%   n'th cell of *c*, `c{m}`, is a 2D matrix of size $M(n) \times W$ and
%   containing the output from the m'th channel subsampled at a rate of
%   $a(m)$.  *c\{m\}(n,l)* is thus the value of the coefficient for time index
%   *n*, frequency index *m* and signal channel *l*.
%
%   The coefficients *c* computed from the signal *f* and the filterbank
%   with windows *g_m* are defined by
%
%   ..            L-1
%      c_m(n+1) = sum f(l+1) * g_m (a(m)n-l+1)
%                 l=0
%
%   .. math:: c_m\left(n+1\right)=\sum_{l=0}^{L-1}f\left(l+1\right)g\left(a_mn-l+1\right)
%
%   where $an-l$ is computed modulo $L$.
%
%   See also: ufilterbank, ifilterbank, pfilt
%
%   References: bohlfe02
    
if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.L=[];
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,'FILTERBANK',0);

mustbeuniform=0;
  
if ~isnumeric(a)
  error('%s: a must be numeric.',upper(callfun));
end;
  
if isempty(L)
  L=filterbanklength(Ls,a);
end;

[g,info]=filterbankwin(g,a,L,'normal');
M=info.M;

if length(a)>1 
  if  length(a)~=M
    error(['%s: The number of entries in "a" must match the number of ' ...
           'filters.'],upper(callfun));
  end;
else
  a=a*ones(M,1);
end;

N=L./a;

c=cell(M,1);
for m=1:M
  c{m}=zeros(N(m),W);
end;
  
G=zeros(L,M);
for ii=1:M
  G(:,ii)=fft(fir2long(g{ii},L));
end;

for w=1:W
  F=fft(f(:,w),L);
  for m=1:M
    c{m}(:,w)=ifft(sum(reshape(F.*G(:,m),N(m),a(m)),2))/a(m);
  end;
end;
