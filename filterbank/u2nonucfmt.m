function c = u2nonucfmt(cu, p)
%U2NONUCFMT Uniform to non-uniform filterbank coefficient format
%   Usage:  c=u2nonucfmt(cu,pk)
%
%   Input parameters:
%         cu   : Uniform filterbank coefficients.
%
%   Output parameters:
%         c    : Non-uniform filterbank coefficients.
%         p    : Numbers of copies of each filter.
%
%   `c = u2nonucfmt(cu,pk)` changes the coefficient format from
%   uniform filterbank coefficients *cu* (`M=sum(p)` channels) to
%   non-uniform coefficients *c* (`numel(p)` channels)  such that each
%   channel of *c* consinst of `p(m)` interleaved channels of *cu*.
%
%   The output *c* is a cell-array in any case.
%
%   See also: nonu2ufilterbank
%
%   References: akkva2003

complainif_notenoughargs(nargin,2,mfilename);

if isempty(cu)
   error('%s: cu must be non-empty.',upper(mfilename));
end

if iscell(cu)
    if any(cellfun(@isempty,cu));
      error('%s: Elements of cu must be non-empty.',upper(mfilename));
    end

    M = numel(cu);
    W = size(cu{1},2);

    Lc = size(cu{1},1);
    if any(Lc~=cellfun(@(cEl)size(cEl,1),cu))
        error('%s: Coefficient subbands do not have an uniform length',...
              upper(mfilename));
    end
elseif isnumeric(cu)
    M = size(cu,2);
    W = size(cu,3);
    Lc = size(cu,1);
else
    error('%s: cu must be a cell array or numeric.',upper(mfilename));
end

if isempty(p) || ~isvector(p)
   error('%s: p must be a non-empty vector.',upper(mfilename));
end

if sum(p) ~= M
    error(['%s: Total number of filters in p does not comply with ',...
           'number of channels'],upper(mfilename));
end

Mnonu = numel(p);
c = cell(Mnonu,1);
p = p(:);
pkcumsum = cumsum([1;p]);
crange = arrayfun(@(pEl,pcEl)pcEl:pcEl+pEl-1,p,pkcumsum(1:end-1),'UniformOutput',0);

if iscell(cu)
    for m=1:Mnonu
        ctmp = [cu{crange{m}}].';
        c{m} = reshape(ctmp(:),W,[]).';
    end
else
    for m=1:Mnonu
        c{m} = zeros(p(m)*Lc,W,assert_classname(cu));
        for w=1:W
            c{m}(:,w) = reshape(cu(:,crange{m},w).',1,[]);
        end
    end
end

% Post check whether there is the same number of coefficients
if sum(cellfun(@(cEl) size(cEl,1),c)) ~= M*Lc
    error(['%s: Invalid number of coefficients in subbands.'],upper(mfilename));
end

