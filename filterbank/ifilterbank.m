function f=ifilterbank(c,g,a,varargin);  
%IFILTERBANK  Filter bank inversion
%   Usage:  f=ifilterbank(c,g,a);
%
%   IFILTERBANK(c,g,a) will synthesize a signal f from the coefficients c
%   using the filters stores in g for a channel subsampling rate of _a (the
%   hop-size). The coefficients has to be in the format returned by
%   either FILTERBANK or UFILTERBANK.
%
%   The filter format for g is the same as for FILTERBANK.
%
%   If perfect reconstruction is desired, the filters must be the duals
%   of the filters used to generate the coefficients. See the help on
%   FILTERBANKDUAL.
%
%   See also: filterbank, ufilterbank, filterbankdual
%
%R  bohlfe02

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.Ls=[];
[flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);

[a,M,longestfilter,lcm_a]=assert_filterbankinput(g,a,1);

if iscell(c)
  Mcoef=numel(c);
else
  Mcoef=size(c,2);
end;

if ~(size(c,2)==Mcoef)
  error(['Mismatch between the size of the input coefficients and the ' ...
         'number of filters.']);
end;

if iscell(c)
  error('Not implemented yet.');
else
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
      f(:,w)=f(:,w)+ifft(repmat(fft(c(:,m,w)),a,1).*conj(G(:,m)));
    end;
  end;
  
end;
  
% Cut or extend f to the correct length, if desired.
if ~isempty(Ls)
  f=postpad(f,Ls);
else
  Ls=L;
end;
  
