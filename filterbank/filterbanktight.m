function gtout=filterbanktight(g,a,varargin);
%FILTERBANKTIGHT  Tight filters
%   Usage:  gt=filterbanktight(g,a);
%
%   `filterabanktight(g,a)` computes the canonical tight filters of *g* for a
%   channel subsampling rate of *a* (hop-size).
%
%   The input and output format of the filters *g* are described in the
%   help of |filterbank|.
%
%   `filterabanktight(g,a,L)` computes canonical tight filters for a system
%   of length *L*. If *L* is not specified, the shortest possible transform
%   length is choosen.
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

[g,info]=filterbankwin(g,a,L,'normal');
M=info.M;

if all(a==a(1))
  % Uniform filterbank, use polyphase representation
  if isempty(L)
      error('%s: You need to specify L.',upper(mfilename));
  end;

  a=a(1);

  % G1 is done this way just so that we can determine the data type.
  G1=comp_transferfunction(g{1},L);
  thisclass=assert_classname(G1);
  G=zeros(L,M,thisclass);
  G(:,1)=G1;
  for ii=2:M
    G(:,ii)=comp_transferfunction(g{ii},L);
  end;
  
  N=L/a;
  
  H=zeros(a,M,thisclass);
  gt=zeros(N,M,thisclass);
  
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
    gtout{m}=cast(gt(:,m),thisclass);
  end;
  
else

  error('Not implemented yet.');  
  
end;
