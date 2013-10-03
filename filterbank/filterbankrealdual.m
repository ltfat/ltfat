function gdout=filterbankrealdual(g,a,varargin);
%FILTERBANKREALDUAL  Dual filters of filterbank for real signals only 
%   Usage:  gd=filterbankdual(g,a);
%
%   `filterabankdual(g,a)` computes the canonical dual filters of *g* for a
%   channel subsampling rate of *a* (hop-size). The dual filters work only
%   for real-valued signals. Use this function on the common construction
%   where the filters in *g* only covers the positive frequencies.
%
%   The format of the filters *g* are described in the
%   help of |filterbank|.
%
%   To actually invert the output of a filterbank, use the dual filters
%   together with `2*real(ifilterbank(...))`.
%
%   See also: filterbank, ufilterbank, ifilterbank

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.L=[];
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

[g,info]=filterbankwin(g,a,L,'normal');
M=info.M;

if (~isempty(L)) && (L~=filterbanklength(L,a))
        error(['%s: Specified length L is incompatible with the length of ' ...
               'the time shifts.'],upper(mfilename));
end;

if info.isuniform
        
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
  
  % This is the original code
  %for k=0:a-1
  %  Ha(k+1,:) =      G(mod(w-k*N,L)+1,:);
  %  Hb(k+1,:) = conj(G(mod(k*N-w,L)+1,:));
  %end;
  
  gd=zeros(N,M,thisclass);
  
  for w=0:N-1
    idx_a = mod(w-(0:a-1)*N,L)+1;
    idx_b = mod((0:a-1)*N-w,L)+1;
    Ha = G(idx_a,:);
    Hb = conj(G(idx_b,:));
    
    Ha=(Ha*Ha'+Hb*Hb')\Ha;
    
    gd(idx_a,:)=Ha;
  end;
  
  gd=ifft(gd)*a;
  
  if isreal(g)
    gd=real(gd);
  end;
  
  gdout=cell(1,M);
  for m=1:M
    gdout{m}=struct('h',cast(gd(:,m),thisclass),'offset',0);
  end;
  
else

    if info.ispainless
        F=comp_filterbankresponse(g,info.a,L,1);
        
        gdout=cell(1,M);
        for m=1:M
            gl=numel(g{m}.H);
            thisgd=struct();
            H=circshift(comp_transferfunction(g{m},L)./F,-g{m}.foff);
            thisgd.H=H(1:gl);
            thisgd.foff=g{m}.foff;
            thisgd.realonly=0;
            thisgd.delay=0;
            
            gdout{m}=thisgd;
        end;
    else
        error(['%s: The canonical dual frame of this system is not a ' ...
               'filterbank. You must call an iterative ' ...
               'method to perform the desired inverstion. Please see ' ...
               'FRANAITER or FRSYNITER.'],upper(mfilename));                
        
    end;  
end;
