function gdout=filterbankdual(g,a,varargin);
%FILTERBANKDUAL  Dual filters
%   Usage:  gd=filterbankdual(g,a);
%
%   FILTERABANKDUAL(g,a) computes the canonical dual filters of g for a
%   channel subsampling rate of _a (hop-size).
%
%   The input and output format of the filters g are described in the
%   help of FILTERANK.
%
%   To actually invert the output of a filterbank, use the dual filters
%   together with the IFILTERBANK function.
%
%   See also: filterbank, ifilterbank

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~iscell(g)
  error('g must be a cell-array.');
end;

M=numel(g);
  
longestfilter=max(cellfun(@numel,g));

definput.keyvals.L=[];
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

if isempty(L)
  L=ceil(longestfilter/a)*a;
end;

G=zeros(L,M);
for ii=1:M
  G(:,ii)=fft(fir2long(g{ii},L));
end;

N=L/a;

H=zeros(a,M);
gd=zeros(N,M);

for w=0:N-1
  for k=0:a-1
    for m=0:M-1
      H(k+1,m+1)=G(mod(w-k*N,L)+1,m+1);
    end;
  end;
      
  H=pinv(H)';

  for k=0:a-1
    for m=0:M-1
      gd(mod(w-k*N,L)+1,m+1)=H(k+1,m+1);
    end;
  end;
end;
  
gd=ifft(gd)*a;

if isreal(g)
  gd=real(gd);
end;

gdout=cell(1,M);
for m=1:M
  gdout{m}=gd(:,m);
end;