function gdout=filterbankdual(g,a,varargin);
%FILTERBANKDUAL  Dual filters
%   Usage:  gd=filterbankdual(g,a);
%           gd=filterbankdual(g,a,L);
%
%   `filterabankdual(g,a)` computes the canonical dual filters of *g* for a
%   channel subsampling rate of *a* (hop-size).
%
%   The input and output format of the filters *g* are described in the
%   help of |filterbank|.
%
%   `filterabankdual(g,a,L)` computes canonical dual filters for a system
%   of length *L*. If *L* is not specified, the shortest possible
%   transform length is choosen.
%
%   To actually invert the output of a filterbank, use the dual filters
%   together with the |ifilterbank| function.
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
else
  if rem(L,lcm_a)>0
    error(['%s: Specified length L is incompatible with the length of ' ...
           'the time shifts. L = %i, lcm_a = %i'],upper(mfilename),L,lcm_a);
  end;
end;


if all(a==a(1))
  % Uniform filterbank, use polyphase representation
  a=a(1);
  
  G=zeros(L,M,assert_classname(g{1}));
  for ii=1:M
    G(:,ii)=fft(fir2long(g{ii},L));
  end;
  
  N=L/a;
  
  gd=zeros(N,M,assert_classname(g{1}));
  
  for w=0:N-1
    idx = mod(w-(0:a-1)*N,L)+1;
    H = G(idx,:);
    
    H=pinv(H)';
    
    gd(idx,:)=H;
  end;
  
  gd=ifft(gd)*a;
  
  if isreal(g)
    gd=real(gd);
  end;
  
  gdout=cell(1,M);
  for m=1:M
    gdout{m}=cast(gd(:,m),assert_classname(g{1}));
  end;
  
else
  
  error('Not implemented yet.');
  
end;
