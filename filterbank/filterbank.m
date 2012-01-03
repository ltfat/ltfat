function c=filterbank(f,g,a);  
%FILTERBANK   Apply filterbank
%   Usage:  c=filterbank(f,g,a);
%
%   FILTERBANK(f,g,a) applies the filters given in g to the signal f. Each
%   subband will be subsampled by a factor of _a (the hop-size). In contrast
%   to UFILTERBANK, _a can be a vector so the hop-size can be
%   channel-dependant. If f is a matrix, the transformation is applied to
%   each column.
%
%   The filters g must be a cell-array, where each entry in the cell
%   array corresponds to an FIR filter.
%
%   The output coefficients are stored a cell array. More precisely, the nth
%   cell of c, c{m}, is a 2D matrix of size M(n) x W and containing the
%   output from the m'th channel subsampled at a rate of a(m).  c{m}(n,l) is
%   thus the value of the coefficient for time index n, frequency index m
%   and signal channel l.
%
%R  bohlfe02
  
if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

[a,M,longestfilter,lcm_a]=assert_filterbankinput(g,a);

[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,'FILTERBANK',0);

L=ceil(max(Ls,longestfilter)/lcm_a)*lcm_a;

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
  F=fft(postpad(f(:,w),L));
  for m=1:M
    c{m}(:,w)=ifft(sum(reshape(F.*G(:,m),N(m),a(m)),2))/a(m);
  end;
end;


%OLDFORMAT
