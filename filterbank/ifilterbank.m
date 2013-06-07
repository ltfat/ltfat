function f=ifilterbank(c,g,a,varargin);  
%IFILTERBANK  Filter bank inversion
%   Usage:  f=ifilterbank(c,g,a);
%
%   `ifilterbank(c,g,a)` synthesizes a signal *f* from the coefficients *c*
%   using the filters stored in *g* for a channel subsampling rate of *a* (the
%   hop-size). The coefficients has to be in the format returned by
%   either |filterbank| or |ufilterbank|.
%
%   The filter format for *g* is the same as for |filterbank|.
%
%   If perfect reconstruction is desired, the filters must be the duals
%   of the filters used to generate the coefficients. See the help on
%   |filterbankdual|.
%
%   See also: filterbank, ufilterbank, filterbankdual
%
%   References: bohlfe02

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.Ls=[];
[flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);

if ~isnumeric(a)
  error('%s: a must be numeric.',upper(callfun));
end;

if iscell(c)
  Mcoef=numel(c);
  W=size(c{1},2);
else
  Mcoef=size(c,2);
  W=size(c,3);    
end;

mustbeuniform=isnumeric(c);

L=filterbanklengthcoef(c,a);

[g,info]=filterbankwin(g,a,L,'normal');

M=info.M;

if length(a)>1 
  if  length(a)~=M
    error(['%s: The number of entries in "a" must match the number of ' ...
           'filters.'],upper(callfun));
  end;
else
  a=a*ones(M,1);
end;

if ~(M==Mcoef)
  error(['Mismatch between the size of the input coefficients and the ' ...
         'number of filters.']);
end;

if iscell(c)
    f=zeros(L,W,assert_classname(c{1}));
else
    a=a(1);
    f=zeros(L,W,assert_classname(c));
end;

l=(0:L-1).'/L;
for m=1:M
    if isfield(g{m},'h')
        g_time=circshift(postpad(g{m}.h,L),g{m}.offset).*...
               exp(2*pi*1i*round(g{m}.centre*L/2)*l);
        conjG=conj(fft(g_time));
    else
        conjG=conj(circshift(postpad(g{m}.H(L),L),g{m}.foff(L)).*exp(-2*pi*1i*round(g{m}.delay)*l));
    end;
    
    if iscell(c)
        for w=1:W                        
            % This repmat cannot be replaced by bsxfun
            f(:,w)=f(:,w)+ifft(bsxfun(@times,repmat(fft(c{m}(:,w)),a(m),1),conjG));
        end;                
    else
        for w=1:W
            % This repmat cannot be replaced by bsxfun
            f(:,w)=f(:,w)+ifft(repmat(fft(c(:,m,w)),a,1).*conjG);
        end;        
    end;
end;

% Old stuff from before the struct filters
if 0
    if iscell(c)
        
        for m=1:M
            for w=1:W
                
                G=fft(fir2long(g{m},L));
                
                % This repmat cannot be replaced by bsxfun
                f(:,w)=f(:,w)+ifft(bsxfun(@times,repmat(fft(c{m}(:,w)),a(m),1),conj(G)));
            end;
        end;
        
    else
        
        a=a(1);
        
        G=zeros(L,M,assert_classname(c,g{1}));
        for ii=1:M
            G(:,ii)=fft(fir2long(g{ii},L));
        end;
        
        f=zeros(L,W,assert_classname(c,g{1}));
        for w=1:W
            for m=1:M
                % This repmat cannot be replaced by bsxfun
                f(:,w)=f(:,w)+ifft(repmat(fft(c(:,m,w)),a,1).*conj(G(:,m)));
            end;
        end;
        
    end;
    
end;
  
% Cut or extend f to the correct length, if desired.
if ~isempty(Ls)
  f=postpad(f,Ls);
else
  Ls=L;
end;
