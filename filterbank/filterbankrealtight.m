function gtout=filterbankrealtight(g,a,varargin);
%FILTERBANKTIGHT  Tight filters of filterbank for real signals only 
%   Usage:  gd=filterbanktight(g,a);
%
%   FILTERABANKTIGHT(g,a) computes the canonical tight filters of g for a
%   channel subsampling rate of _a (hop-size). The tight filters work only
%   for real-valued signals. Use this function on the common construction
%   where the filters in g only covers the positive frequencies.
%
%   The format of the filters g are described in the
%   help of FILTERANK.
%
%   To actually invert the output of a filterbank, use the tight filters
%   together with 2*real(IFILTERBANK(...)) function.
%
%   See also: filterbank, ufilterbank, ifilterbank

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

  gt=zeros(N,M);
  
  for w=0:N-1
    idx_a = mod(w-(0:a-1)*N,L)+1;
    idx_b = mod((0:a-1)*N-w,L)+1;
    Ha = G(idx_a,:);
    Hb = conj(G(idx_b,:));
    
    Ha=sqrtm(Ha*Ha'+Hb*Hb')\Ha;
    
    gt(idx_a,:)=Ha;
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



