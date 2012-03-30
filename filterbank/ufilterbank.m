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
%   See also: ifilterbank, filterbankdual
%
%   References: bohlfe02
  
if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.L=[];
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,'UFILTERBANK',0);

a=a(1);

if isempty(L)
  L=filterbanklengthsignal(Ls,a);
end;

[g,info]=filterbankwin(g,a,L,'normal');
M=info.M;

N=L/a;

f=postpad(f,L);

gw=zeros(L,M);
for ii=1:M
  gw(:,ii)=fir2long(g{ii},L);
end;

c=comp_ufilterbank_fft(f,gw,a);
