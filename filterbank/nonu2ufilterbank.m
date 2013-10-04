function [gu,au,info]=nonu2ufilterbank(g,a)
%NONU2UFILTERBANK   Non-uniform to Uniform filterbank transform
%   Usage:  [gu,au]=nonu2ufilterbank(g,a)
%
%   `[gu,au]=nonu2ufilterbank(g,a)` calculates uniform filterbank *gu*, 
%   `au=lcm(a)` which is identical to the (possibly non-uniform) filterbank
%   *g*,*a* in terms of the equal output coefficients. Each filter `g{k}` 
%   is replaced by $p=au/a(k)$ delayed versions of itself
%   $z^{-ma(k)}G_k(z)$ for $m=0,\ldots,p-1$.
%
%   This allows using the factorisation algorithm when determining
%   filterbank frame bounds in |filterbankbounds| and
%   |filterbankrealbounds| and in the computation of the dual filterbank 
%   in |filterbankdual| and |filterbakrealdual| which do not work 
%   with non-uniform filterbanks.
%
%   See also: ufilterbank, filterbank, filterbankbounds, filterbankdual
%

%   References: TBD

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(a)
  error('%s: a must be numeric.',upper(mfilename));
end;

if ~iscell(g) || ...
   ~all(cellfun(@(gEl) isstruct(gEl) && (isfield(gEl,'h')||isfield(gEl,'H')),g))
  error('%s: a must be a cell array of structs containing filter definition.',upper(mfilename));
end;

a = comp_filterbank_a(a,numel(g));

if size(a,2)==2 && ~all(a(:,2)==1) && rem(a(:,1),1)~=0
   error('%s: Filterbanks with fractional subsampling are not supported.',upper(mfilename)); 
end

au=filterbanklength(1,a);

pk = au./a(:,1);

gu=cell(sum(pk),1);

auIdx = 1;
for m=1:numel(g)
   for ii=0:pk(m)-1
      gu{auIdx} = g{m};
      if(isfield(gu{auIdx},'H'))
         if(~isfield(gu{auIdx},'delay'))
            gu{auIdx}.delay = 0;
         end
         gu{auIdx}.delay = gu{auIdx}.delay-a(m)*ii;
      end
      
      if(isfield(gu{auIdx},'h'))
         if(~isfield(gu{auIdx},'offset'))
            gu{auIdx}.offset = 0;
         end
         gu{auIdx}.offset = gu{auIdx}.offset-a(m)*ii;
      end
      auIdx = auIdx+1;
   end
end

if nargout>2
   % Array of lengths of identical filterbanks of each channel
   info.p = pk;
end
