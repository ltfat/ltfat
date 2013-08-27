function c=ufilterbank(f,g,a,varargin);  
%UFILTERBANK   Apply Uniform filterbank
%   Usage:  c=ufilterbank(f,g,a);
%
%   `ufilterbank(f,g,a)` applies the filter given in *g* to the signal
%   *f*. Each subband will be subsampled by a factor of *a* (the
%   hop-size). If *f* is a matrix, the transformation is applied to each
%   column.
%
%   The filters *g* must be a cell-array, where each entry in the cell
%   array corresponds to an FIR filter.
%
%   If *f* is a single vector, then the output will be a matrix, where each
%   column in *f* is filtered by the corresponding filter in *g*. If *f* is
%   a matrix, the output will be 3-dimensional, and the third dimension will
%   correspond to the columns of the input signal.
%
%   The coefficients *c* computed from the signal *f* and the filterbank
%   with windows *g_m* are defined by
%
%   ..              L-1
%      c(n+1,m+1) = sum f(l+1) * g_m (an-l+1)
%                   l=0
%
%   .. math:: c\left(n+1,m+1\right)=\sum_{l=0}^{L-1}f\left(l+1\right)g\left(an-l+1\right)

%
%   See also: ifilterbank, filterbankdual
%
%   References: bohlfe02
  
if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'pfilt'};
definput.keyvals.L=[];
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,'UFILTERBANK',0);

a=a(1);

if isempty(L)
  L=filterbanklength(Ls,a);
end;

[g,info]=filterbankwin(g,a,L,'normal');
M=info.M;
N=L/a;

f=postpad(f,L);

c=zeros(N,M,W,assert_classname(f));
g = comp_filterbank_pre(g,info.a,L,kv.crossover);

ctmp=comp_filterbank(f,g,info.a);


for m=1:M    
    c(:,m,:)=ctmp{m};
end;
