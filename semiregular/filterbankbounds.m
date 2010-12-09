function [AF,BF]=filterbankbounds(g,a,varargin);
%FILTERBANKBOUNDS  Frame bounds of filter bank
%   Usage: fcond=filterbankbounds(g,a);
%          [A,B]=filterbankbounds(g,a);
%
%   FILTERBANKBOUNDS(g,a) calculates the ratio B/A of the frame bounds of
%   the filterbank specified by g and _a. The ratio is a measure of the
%   stability of the system.
%
%   [A,B]=FILTERBANKBOUNDS(g,a) returns the lower and upper frame bounds
%   explicitly.
%
%   See also: filterbank, filterbankdual
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

[M,longestfilter]=assert_filterbankinput(g,a);

definput.keyvals.L=[];
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

if isempty(L)
  L=ceil(longestfilter/a)*a;
end;

N=L/a;

G=zeros(L,M);
for ii=1:M
  G(:,ii)=fft(fir2long(g{ii},L));
end;

H=zeros(a,M);

AF=Inf;
BF=0;

for w=0:N-1
  
  for k=0:a-1
    for m=0:M-1
      H(k+1,m+1)=G(mod(w-k*N,L)+1,m+1);
    end;
  end;
  
  % A 'real' is needed here, because the matrices are known to be
  % Hermitian, but sometimes Matlab/Octave does not recognize this.
  work=real(eig(H*H'));
  AF=min(AF,min(work));
  BF=max(BF,max(work));
  
end;
    
AF=AF/a;
BF=BF/a;

if nargout<2
  % Avoid the potential warning about division by zero.
  if AF==0
    AF=Inf;
  else
    AF=BF/AF;
  end;
end;
  