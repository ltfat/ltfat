function c=pfilt(f,g,varargin)
%PFILT  Apply filter with periodic boundary conditions
%   Usage:  h=pfilt(f,g);
%           h=pfilt(f,g,a,dim);
%
%   `pfilt(f,g)` applies the filter *g* to the input *f*. If *f* is a
%   matrix, the filter is applied along each column.
%
%   `pfilt(f,g,a)` does the same, but downsamples the output keeping only
%   every a'th sample (starting with the first one).
%
%   `pfilt(f,g,a,dim)` filters along dimension dim. The default value of
%   [] means to filter along the first non-singleton dimension.
%
%   The filter *g* can be a vector, in which case the vector is treated
%   as a zero-delay FIR filter.
%
%   The filter *g* can be a cell array. The following options are
%   possible:
%
%     * If the first element of the cell array is the name of one of the
%       windows from |firwin|, the whole cell array is passed onto
%       |firfilter|.
%
%     * If the first element of the cell array is `'bl'`, the rest of the
%       cell array is passed onto |blfilter|.
%
%     * If the first element of the cell array is `'pgauss'`, `'psech'`,
%       the rest of the parameters is passed onto the respective
%       function. Note that you do not need to specify the length *L*.
%
%   The coefficients obtained from filtering a signal *f* by a filter *g* are
%   defined by
%
%   ..          L-1
%      c(n+1) = sum f(l+1) * g(an-l+1)
%               l=0
%
%   .. math:: c\left(n+1\right)=\sum_{l=0}^{L-1}f\left(l+1\right)g\left(an-l+1\right)
%
%   where $an-l$ is computed modulo $L$.
%
%   See also: pconv

  
% Assert correct input.
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'pfilt'};
definput.keyvals.a=1;
definput.keyvals.dim=[];
[flags,kv,a,dim]=ltfatarghelper({'a','dim'},definput,varargin);

[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,[],dim,upper(mfilename));

[g,info] = comp_fourierwindow(g,L,upper(mfilename));

outIsReal = isreal(f) &&...
              (isfield(g,'h') && isreal(g.h) && ~(isfield(g,'fc') && g.fc~=0) || isfield(g,'H') && g.realonly);

asan = comp_filterbank_a(a,1);
g = comp_filterbank_pre({g},a,L,kv.crossover);

c = comp_filterbank(f,g,a);
c = c{1};

permutedsize(1)=size(c,1);
  
c=assert_sigreshape_post(c,dim,permutedsize,order);

if outIsReal
   c = real(c);
end

