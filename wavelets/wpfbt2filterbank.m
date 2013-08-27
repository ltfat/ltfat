function [g,a] = wpfbt2filterbank( wtdef, varargin)
%WPFBT2FILTERBANK  WPFBT equivalent non-iterated filterbank
%   Usage: [g,a] = wpfbt2filterbank(wtdef)
%
%   Input parameters:
%         wtdef : Wavelet filter tree definition
%
%   Output parameters:
%         g   : Cell array containing filters
%         a   : Vector of sub-/upsampling factors
%
%   `wpfbtmultid(wtdef)` calculates the impulse responses *g* and the 
%   subsampling factors *a* of non-iterated filterbank, which is equivalent
%   to the wavelet filterbank tree described by *wtdef*. The returned 
%   parameters can be used directly in |filterbank|, |ufilterbank| or 
%   |filterbank|. 
%   
%   The filters are scaled if *a* is not returned. 
%   
%   See also: wfbtinit, wfbt2filterbank, filterbank


if(nargin<1)
    error('%s: Not enough input arguments',upper(mfilename));
end

% build the tree
wt = wfbtinit({'strict',wtdef},varargin{:});

wtPath = nodesBForder(wt);
rangeLoc = cellfun(@(nEl) 1:numel(nEl.g),wt.nodes(wtPath),'UniformOutput',0);
rangeOut = mat2cell(1:sum(rangeLocCount),1,...
           cellfun(@(nEl) numel(nEl.g),wt.nodes(wtPath)));


[g,a] = nodesMultid(wtPath,rangeLoc,rangeOut,wt);

if nargout<2
   % Scale filters if a is not returned
   for nn=1:numel(g)
       g{nn}.h = g{nn}.h/sqrt(a(nn));
   end
end

