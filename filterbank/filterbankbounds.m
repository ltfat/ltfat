function [AF,BF]=filterbankbounds(g,a,varargin);

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~iscell(g)
  error('g must be a cell-array.');
end;

definput.keyvals.L=[];
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

M=numel(g);

longestfilter=max(cellfun(@numel,g));

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