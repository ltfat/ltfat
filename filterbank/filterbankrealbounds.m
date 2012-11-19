function [AF,BF]=filterbankrealbounds(g,a,varargin);
%FILTERBANKREALBOUNDS  Frame bounds of filter bank for real signals only
%   Usage: fcond=filterbankrealbounds(g,a);
%          [A,B]=filterbankrealbounds(g,a);
%
%   `filterbankrealbounds(g,a)` calculates the ratio $B/A$ of the frame
%   bounds of the filterbank specified by *g* and *a*. The ratio is a measure
%   of the stability of the system.  Use this function on the common
%   construction where the filters in *g* only covers the positive
%   frequencies.
%
%   `[A,B]=filterbankrealbounds(g,a)` returns the lower and upper frame
%   bounds explicitly.
%
%   See also: filterbank, filterbankdual
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

[a,M,longestfilter,lcm_a]=assert_filterbankinput(g,a);


definput.keyvals.L=[];
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

if isempty(L)
  L=ceil(longestfilter/lcm_a)*lcm_a;
end;

AF=Inf;
BF=0;
  
if all(a==a(1))
  % Uniform filterbank, use polyphase representation
  a=a(1);
  
  N=L/a;

  G=zeros(L,M);
  for ii=1:M
    G(:,ii)=fft(fir2long(g{ii},L));
  end;
  
  Ha=zeros(a,M);
  Hb=zeros(a,M);
  
  for w=0:N-1
    idx_a = mod(w-(0:a-1)*N,L)+1;
    idx_b = mod((0:a-1)*N-w,L)+1;
    Ha = G(idx_a,:);
    Hb = conj(G(idx_b,:));
    
    % A 'real' is needed here, because the matrices are known to be
    % Hermitian, but sometimes Matlab/Octave does not recognize this.  
    work=real(eig(real(Ha*Ha'+Hb*Hb')));
    
    AF=min(AF,min(work));
    BF=max(BF,max(work));
    
  end;
  
  AF=AF/a;
  BF=BF/a;
  
else

    error('Not implemented yet.');
  
end;

if nargout<2
  % Avoid the potential warning about division by zero.
  if AF==0
    AF=Inf;
  else
    AF=BF/AF;
  end;
end;
  

