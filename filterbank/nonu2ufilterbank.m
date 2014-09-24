function [gu,au,p]=nonu2ufilterbank(g,a)
%NONU2UFILTERBANK   Non-uniform to uniform filterbank transform
%   Usage:  [gu,au]=nonu2ufilterbank(g,a)
%
%   Input parameters:
%         g     : Filters as a cell array of structs.
%         a     : Subsampling factors.
%
%   Output parameters:
%         gu    : Filters as a cell array of structs.
%         au    : Uniform subsampling factor.
%         pk    : Numbers of copies of each filter.
%
%   `[gu,au]=nonu2ufilterbank(g,a)` calculates uniform filterbank `gu`, 
%   `au=lcm(a)` which is identical to the (possibly non-uniform) filterbank
%   *g*, *a* in terms of the equal output coefficients. Each filter `g{k}` 
%   is replaced by $p(k)=au/a(k)$ advanced versions of itself such that
%   $z^{ma(k)}G_k(z)$ for $m=0,\ldots,p-1$.
%
%   This allows using the factorisation algorithm when determining
%   filterbank frame bounds in |filterbankbounds| and
%   |filterbankrealbounds| and in the computation of the dual filterbank 
%   in |filterbankdual| and |filterbankrealdual| which do not work 
%   with non-uniform filterbanks.
%
%   One can change between the coefficient formats of `gu`, `au` and 
%   `g`, `a` using |nonu2ucfmt| and |u2nonucfmt| in the reverse direction.
%
%   See also: ufilterbank, filterbank, filterbankbounds, filterbankdual
%
%   References: akkva2003

complainif_notenoughargs(nargin,2,'NONU2UFILTERBANK');

if ~isnumeric(a) || isempty(a)
  error('%s: a must be numeric.',upper(mfilename));
end;

% Allow only structs and numeric arrays as filters
if isempty(g) || ~iscell(g) || any(cellfun(@isempty,g)) || ...
   ~all(cellfun(@(gEl) isnumeric(gEl) || ...
   (isstruct(gEl) && (isfield(gEl,'h')||isfield(gEl,'H'))),g))
  error(['%s: g must be a cell array of structs containing filter',...
         ' definition or numeric arrays.'],upper(mfilename));
end;

% Do not allow *g* in format {'dual',g} etc. since it requires L and
% we want to work without it here.
% This case is also captured by the next check, but we want a separate 
% error message.
if ischar(g{1})
   error('%s: %s option is not supported in this function.',...
         upper(mfilename),g{1}); 
end

% Sanitize a
a = comp_filterbank_a(a,numel(g));

% This does not work for fractional subsampling
if size(a,2)==2 && ~all(a(:,2)==1) && rem(a(:,1),1)~=0
   error(['%s: Filterbanks with fractional subsampling are not',...
          ' supported.'],upper(mfilename)); 
end

% This does lcm(a)
au=filterbanklength(1,a);

% Numbers of copies of each filter
p = au./a(:,1);

if all(a(:,1)==a(1,1))
    % Filterbank is already uniform, there is nothing to be done.
    gu = g;
    au = a;
    return;
end

% Sanitize filters which come as numeric vertors.
gIdx = cellfun(@isnumeric,g);
if any(gIdx)
    g(gIdx) = filterbankwin(g(gIdx),a,'normal');
end

% Do the actual filter copies
gu=cell(sum(p),1);
auIdx = 1;
for m=1:numel(g)
   for ii=0:p(m)-1
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


