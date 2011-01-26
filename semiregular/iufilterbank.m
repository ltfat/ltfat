function f=iufilterbank(c,g,a);  
%IUFILTERBANK  Filter bank inversion
%   Usage:  f=iufilterbank(c,g,a);
%
%   IUFILTERBANK(c,g,a) will synthesize a signal f from the coefficients c
%   using the filters stores in g for a channel subsampling rate of _a (the
%   hop-size).
%
%   The filter format for g is the same as for FILTERBANK.
%
%   If perfect reconstruction is desired, the filters must be the duals
%   of the filters used to generate the coefficients. See the help on
%   FILTERBANKDUAL.
%
%   See also: filterbank, filterbankdual
%
%R  bohlfe02

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

[M,longestfilter]=assert_filterbankinput(g,a);

if ~(size(c,2)==M)
  error(['Mismatch between the size of the input coefficients and the ' ...
         'number of filters.']);
end;

N=size(c,1);
L=N*a;
W=size(c,3);

G=zeros(L,M);
for ii=1:M
  G(:,ii)=fft(fir2long(g{ii},L));
end;

f=zeros(L,W);
for w=1:W
  for m=1:M
    f(:,w)=f(:,w)+ifft(repmat(fft(c(:,m,w)),a,1).*G(:,m));
  end;
end;


  
