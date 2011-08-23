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
%   See also: filterbank, filterbankdual, ufilterbank, ifilterbank

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

[a,M,longestfilter,lcm_a]=assert_filterbankinput(g,a);

definput.keyvals.L=[];
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

if isempty(L)
  L=ceil(longestfilter/lcm_a)*lcm_a;
end;

if all(a==a(1))
  % Uniform filterbank, use polyphase representation
  a=a(1);
  
  G=zeros(L,M);
  for ii=1:M
    G(:,ii)=fft(fir2long(g{ii},L));
  end;
  
  N=L/a;
  
  H=zeros(a,M);
  gt=zeros(N,M);
  
  for w=0:N-1
    idx = mod(w-(0:a-1)*N,L)+1;
    H = G(idx,:);
    
    [U,S,V]=svd(H,'econ');
    H=U*V';  
    
    gt(idx,:)=H;
  end;
  
  gt=ifft(gt)*sqrt(a);  

  if isreal(g)
    gt=real(gt);
  end;
  
  gtout=cell(1,M);
  for m=1:M
    gtout{m}=gt(:,m);
  end;
  
else

  error('Not implemented yet.');  
  
end;