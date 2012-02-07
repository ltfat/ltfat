function h=pfilt(f,g,varargin)
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
%   See also: pconv

  
% Assert correct input.
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.a=1;
definput.keyvals.dim=[];
[flags,kv,a,dim]=ltfatarghelper({'a','dim'},definput,varargin);

L=[];

[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,L,dim,'PFILT');

[g,info] = comp_fourierwindow(g,L,'PFILT');

h=squeeze(comp_ufilterbank_fft(f,g,a));

% FIXME: This check should be removed when comp_ufilterbank_fft{.c/.cc}
% have been fixed.
if isreal(f) && isreal(g)
  h=real(h);
end;

permutedsize(1)=size(h,1);
  
h=assert_sigreshape_post(h,dim,permutedsize,order);

