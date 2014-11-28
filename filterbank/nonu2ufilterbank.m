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

try
    [g,asan] = filterbankwin(g,a,'normal');
catch
    err = lasterror;
    if strcmp(err.identifier,'L:undefined')
        % If it blotched because of the undefined L, explain that.
        % This should capture only formats like {'dual',...} and {'gauss'}
        error(['%s: Function cannot handle g in a format which ',...
               'requires L. Consider pre-formatting the filterbank by ',...
               'calling g = FILTERBANKWIN(g,a,L).'],upper(mfilename));
    else
        % Otherwise just rethrow the error
        error(err.message);
    end
end

% This function does not work for fractional subsampling
if size(asan,2)==2 && ~all(asan(:,2)==1) && rem(asan(:,1),1)~=0
   error(['%s: Filterbanks with fractional subsampling are not',...
          ' supported.'],upper(mfilename)); 
end

% This is effectivelly lcm(a)
au=filterbanklength(1,asan);

% Numbers of copies of each filter
p = au./asan(:,1);

if all(asan(:,1)==asan(1,1))
    % Filterbank is already uniform, there is nothing to be done.
    gu = g;
    au = asan;
    return;
end

% Do the actual filter copies
% This only changes .delay or .offset
gu=cell(sum(p),1);
auIdx = 1;
for m=1:numel(g)
   for ii=0:p(m)-1
      gu{auIdx} = g{m};
      if(isfield(gu{auIdx},'H'))
         if(~isfield(gu{auIdx},'delay'))
            gu{auIdx}.delay = 0;
         end
         gu{auIdx}.delay = gu{auIdx}.delay-asan(m)*ii;
      end
      
      if(isfield(gu{auIdx},'h'))
         if(~isfield(gu{auIdx},'offset'))
            gu{auIdx}.offset = 0;
         end
         gu{auIdx}.offset = gu{auIdx}.offset-asan(m)*ii;
      end
      auIdx = auIdx+1;
   end
end


