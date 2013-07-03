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

L=filterbanklengthcoef(c,a);

[g,info]=filterbankwin(g,a,L,'normal');
M=info.M;
a=info.a;

if iscell(c)
  Mcoef=numel(c);
  W=size(c{1},2);
else
  Mcoef=size(c,2);
  W=size(c,3);    
end;

if ~(M==Mcoef)
  error(['Mismatch between the size of the input coefficients and the ' ...
         'number of filters.']);
end;

if iscell(c)
    f=zeros(L,W,assert_classname(c{1}));
else
    f=zeros(L,W,assert_classname(c));
end;

% This routine must handle the following cases
%
%   * Time-side or frequency-side filters (test for  isfield(g,'H'))
%
%   * Cell array or matrix input (test for iscell(c))
%
%   * Regular or fractional subsampling (test for info.isfractional)

l=(0:L-1).'/L;
for m=1:M
    conjG=conj(comp_transferfunction(g{m},L));
        
    % Handle fractional subsampling (this implies frequency side filters)
    if info.isfractional
        N=size(c{m},1);
        Llarge=ceil(L/N)*N;
        amod=Llarge/N;
        
        for w=1:W                        
            % This repmat cannot be replaced by bsxfun
            innerstuff=middlepad(circshift(repmat(fft(c{m}(:,w)),amod,1),-g{m}.foff),L);
            f(:,w)=f(:,w)+ifft(circshift(bsxfun(@times,innerstuff,circshift(conjG,-g{m}.foff)),g{m}.foff));
        end;                
    else
        if iscell(c)
            for w=1:W
                % This repmat cannot be replaced by bsxfun
                f(:,w)=f(:,w)+ifft(repmat(fft(c{m}(:,w)),a(m),1).*conjG);
            end;
        else
            for w=1:W
                % This repmat cannot be replaced by bsxfun
                f(:,w)=f(:,w)+ifft(repmat(fft(c(:,m,w)),a(m),1).*conjG);
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
