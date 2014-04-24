function cu = nonu2ucfmt(c, p)
%NONU2UCFMT Non-uniform to uniform filterbank coefficient format
%   Usage:  cu=nonu2ucfmt(c,pk)
%
%   Input parameters:
%         c   : Non-uniform filterbank coefficients.
%
%   Output parameters:
%         cu  : Uniform filterbank coefficients.
%         p   : Numbers of copies of each filter.
%
%   `cu = nonu2ucfmt(c,p)` changes the coefficient format from
%   non-uniform filterbank coefficients *c* (`M=numel(p)` channels) to
%   uniform coefficients *c* (`sum(p)` channels)  such that each
%   channel of *cu* consinst of de-interleaved samples of channels of *c*.
%
%   The output *cu* is a cell-array in any case.
%
%   See also: nonu2ufilterbank, u2nonucfmt
%
%   References: akkva2003

complainif_notenoughargs(nargin,2,mfilename);

if isempty(c)
   error('%s: c must be non-empty.',upper(mfilename));
end

if isempty(p) || ~isvector(p)
   error('%s: pk must be a non-empty vector.',upper(mfilename));
end

if iscell(c)
    M = numel(c);
    Lc = cellfun(@(cEl) size(cEl,1),c);
    if all(Lc==Lc(1))
        if ~all(p==1) || numel(p)~=M
           error('%s: Bad format of p for uniform coefficients.',...
              upper(mfilename));
        else
            cu = c;
            % End here, this is already uniform.
            return;
        end
    end
elseif isnumeric(c)
    M = size(c,2);
    if ~all(p==1) || numel(p)~=M
        error('%s: Bad format of p for uniform coefficients.',...
               upper(mfilename));
    end
    % Just convert to cell-array and finish
    cu = cell(M,1);
    for m=1:M
        cu{m}=squeeze(c(:,m,:));
    end;
    % End here, there is nothing else to do.
    return;
else
    error('%s: c must be a cell array or numeric.',upper(mfilename));
end

if numel(p) ~= M
    error(['%s: Number of elements of p does not comply with ',...
           'number of channels passed.'],upper(mfilename));
end


p = p(:);
Mu = sum(p);
cu = cell(Mu,1);

pkcumsum = cumsum([1;p]);
crange = arrayfun(@(pEl,pcEl)pcEl:pcEl+pEl-1,p,pkcumsum(1:end-1),...
                  'UniformOutput',0);

% c can be only cell array at this point
for m=1:M
   for k=1:p(m)
      cu{crange{m}(k)} = c{m}(k:p(m):end,:);
   end
end

% Post check whether the output is really uniform and the numbers of
% coefficients are equal
Lcu = cellfun(@(cEl) size(cEl,1),cu);
if any(Lcu~=Lcu(1)) || sum(Lcu)~=sum(Lc)
    error(['%s: The combination of c and p does not result in uniform ',...
           'coefficients.'],upper(mfilename));
end


