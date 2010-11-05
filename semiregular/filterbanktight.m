function gtout=filterbanktight(g,a,varargin);
%FILTERBANKTIGHT  Tight filters
%   Usage:  gt=filterbanktight(g,a);
%
%   FILTERABANKTIGHT(g,a) computes the canonical tight filters of g for a
%   channel subsampling rate of _a (hop-size).
%
%   The input and output format of the filters g are described in the
%   help of FILTERANK.
%
%   See also: filterbank, ifilterbank, filterbankdual

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
gt=zeros(N,M);

for w=0:N-1
  for k=0:a-1
    for m=0:M-1
      H(k+1,m+1)=G(mod(w-k*N,L)+1,m+1);
    end;
  end;
      
  [U,S,V]=svd(H,0);
  H=U*V';  

  for k=0:a-1
    for m=0:M-1
      gt(mod(w-k*N,L)+1,m+1)=H(k+1,m+1);
    end;
  end;
end;
  
gt=ifft(gt)*sqrt(a);  

if isreal(g)
  gt=real(gt);
end;

gtout=cell(1,M);
for m=1:M
  gtout{m}=gt(:,m);
end;