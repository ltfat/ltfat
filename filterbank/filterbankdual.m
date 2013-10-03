function gdout=filterbankdual(g,a,L,varargin);
%FILTERBANKDUAL  Dual filters
%   Usage:  gd=filterbankdual(g,a);
%           gd=filterbankdual(g,a,L);
%
%   `filterbankdual(g,a)` computes the canonical dual filters of *g* for a
%   channel subsampling rate of *a* (hop-size).
%
%   The input and output format of the filters *g* are described in the
%   help of |filterbank|.
%
%   `filterbankdual(g,a,L)` computes canonical dual filters for a system
%   of length *L*. If *L* is not specified, the shortest possible
%   transform length is choosen.
%
%   To actually invert the output of a filterbank, use the dual filters
%   together with the |ifilterbank| function.
%
%   See also: filterbank, ufilterbank, ifilterbank

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

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
  
  gd=zeros(N,M,thisclass);
  
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
    gdout{m}=struct('h',cast(gd(:,m),thisclass),'offset',0);
  end;
  
else

    if info.ispainless
                
        if isempty(L)
            error('%s: You need to specify L.',upper(mfilename));
        end;

        F=comp_filterbankresponse(g,info.a,L,0);
        
        gdout=cell(1,M);
        for m=1:M
            thisgd=struct();
            H=circshift(comp_transferfunction(g{m},L)./F,-g{m}.foff);
            thisgd.H=H(1:numel(g{m}.H));
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
