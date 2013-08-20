function gtout=filterbankrealtight(g,a,varargin);
%FILTERBANKREALTIGHT  Tight filters of filterbank for real signals only 
%   Usage:  gd=filterbankrealtight(g,a);
%
%   `filterabanktight(g,a)` computes the canonical tight filters of *g* for a
%   channel subsampling rate of *a* (hop-size). The tight filters work only
%   for real-valued signals. Use this function on the common construction
%   where the filters in *g* only covers the positive frequencies.
%
%   The format of the filters *g* are described in the
%   help of |filterbank|.
%
%   To actually invert the output of a filterbank, use the tight filters
%   together with the `ifilterbank` function as in
%   `2*real(ifilterbank(...))`.
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

  gt=zeros(N,M,thisclass);
  
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
    gtout{m}=cast(gt(:,m),thisclass);
  end;
  
else
        
    if info.ispainless
                
        Fsqrt=sqrt(comp_filterbankresponse(g,info.a,L,1));
        
        gtout=cell(1,M);
        for m=1:M
            gl=numel(g{m}.H);
            thisgt=struct();
            H=circshift(comp_transferfunction(g{m},L)./Fsqrt,-g{m}.foff);
            thisgt.H=H(1:gl);
            thisgt.foff=g{m}.foff;
            thisgt.realonly=0;
            thisgt.delay=0;
            
            gtout{m}=thisgt;
        end;
        
    else
        error(['%s: The canonical dual frame of this system is not a ' ...
               'filterbank. You must call an iterative ' ...
               'method to perform the desired inverstion. Please see ' ...
               'FRANAITER or FRSYNITER.'],upper(mfilename));        

    end;
  
end;
