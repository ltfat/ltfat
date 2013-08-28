function [g,a] = wfbt2filterbank( wtdef, varargin)
%WFBT2FILTERBANK  WFBT equivalent non-iterated filterbank
%   Usage: [g,a] = wfbt2filterbank(wtdef)
%
%   Input parameters:
%         wtdef : Wavelet filter tree definition
%
%   Output parameters:
%         g   : Cell array containing filters
%         a   : Vector of sub-/upsampling factors
%
%   `wfbtmultid(wtdef)` calculates the impulse responses *g* and the 
%   subsampling factors *a* of non-iterated filterbank, which is equivalent
%   to the wavelet filterbank tree described by *wtdef*. The returned 
%   parameters can be used directly in |filterbank|, |ufilterbank| or 
%   |filterbank|. 
%   
%   The filters are scaled if *a* is not returned. 
%
%   The function internally calls |wfbtinit| and passes *wtdef* and all 
%   additional parameters to it.   
%   
%   Examples:
%   --------- 
%   
%   The following two examples create a multirate identity filterbank
%   using a tree of depth 3. In the first example, the filterbank is
%   identical to the DWT tree:::
%
%     [g,a] = wfbt2filterbank({'db10',3,'dwt'});
%     
%
%   In the second example, the filterbank is identical to the full
%   wavelet tree:::
%
%     [g,a] = wfbt2filterbank({'db10',3,'full'});
%
%   See also: wfbtinit


if(nargin<1)
    error('%s: Not enough input arguments',upper(mfilename));
end

% build the tree
wt = wfbtinit({'strict',wtdef},varargin{:});

% Pick just nodes with outputs
wtPath = 1:numel(wt.nodes);
wtPath(noOfNodeOutputs(1:numel(wt.nodes),wt)==0)=[];


rangeLoc = rangeInLocalOutputs(wtPath,wt);
rangeOut = rangeInOutputs(wtPath,wt); 

[g,a] = nodesMultid(wtPath,rangeLoc,rangeOut,wt);

if nargout<2
   % Scale filters if a is not returned
   for nn=1:numel(g)
       g{nn}.h = g{nn}.h/sqrt(a(nn));
   end
end
















