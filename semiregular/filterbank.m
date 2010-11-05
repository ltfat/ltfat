function c=filterbank(f,g,a);  
%FILTERBANK   Apply filterbank
%   Usage:  c=filterbank(f,g,a);
%
%   FILTERBANK(f,g,a) applies the filter given in g to the signal f. Each
%   subband will be subsampled by a factor of _a (the hop-size). If f is a
%   matrix, the transformation is applied to each column.
%
%   The filters g must be a cell-array, where each entry in the cell
%   array corresponds to an FIR filter.
%
%   If f is a single vector, then the output will be a matrix, where each
%   column is f filtered by the corresponding filter in g. If f is a
%   matrix, the output will be 3-dimensional, and the third dimension
%   will correspond to the columns of the input signal.
%
%   See also: ifilterbank, filterbankdual
  
if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~iscell(g)
  error('g must be a cell-array.');
end;

M=numel(g);

longestfilter=max(cellfun(@numel,g));

[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,'FILTERBANK',0);

L=ceil(max(Ls,longestfilter)/a)*a;

N=L/a;

c=zeros(N,M,W);
  
G=zeros(L,M);
for ii=1:M
  G(:,ii)=fft(fir2long(g{ii},L));
end;

for w=1:W
  F=fft(f(:,w));
  for m=1:M
    c(:,m,w)=ifft(sum(reshape(F.*conj(G(:,m)),N,a),2))/a;
  end;
end;

  