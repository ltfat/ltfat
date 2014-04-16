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
%   `wpfbt2filterbank(wtdef)` calculates the impulse responses *g* and the 
%   subsampling factors *a* of non-iterated filterbank, which is equivalent
%   to the wavelet packet filterbank tree described by *wtdef*. The returned 
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
%     [g,a] = wpfbt2filterbank({'db10',3,'dwt'});
%     filterbankresponse(g,a,1024,'plot','individual');
%     
%
%   In the second example, the filterbank is identical to the full
%   wavelet tree:::
%
%     [g,a] = wpfbt2filterbank({'db10',3,'full'});
%     filterbankresponse(g,a,1024,'plot','individual');
%
%   See also: wfbtinit


complain_notenoughargs(nargin,1,'WPFBT2FILTERBANK');

definput.import = {'wfbtcommon'};
definput.flags.scale = {'scale','noscale'};
[flags]=ltfatarghelper({},definput,varargin);

% build the tree
wt = wfbtinit({'strict',wtdef},flags.forder);

wtPath = nodesBForder(wt);
rangeLoc = rangeInLocalOutputs(wtPath,wt);

if flags.do_scale
    wtNodes = wt.nodes(wtPath);

    % Scale filters
    for ii=1:numel(wtNodes)
        range = 1:numel(wtNodes{ii}.h);
        range(rangeLoc{ii}) = [];
        wtNodes{ii}.g(range) = ...
            cellfun(@(hEl) setfield(hEl,'h',hEl.h/sqrt(2)),wtNodes{ii}.g(range),...
            'UniformOutput',0);
    end
    % Alter back the tree
    wt.nodes(wtPath) = wtNodes;
end

nIdx = nodesLevelsBForder(wt);
% Now we need to walk the tree by levels
g = {};
a = [];
for ii=1:numel(nIdx)
    rangeLoc = cellfun(@(eEl) 1:numel(eEl.h),wt.nodes(nIdx{ii}),'UniformOutput',0);
    rangeOut = cellfun(@(eEl) numel(eEl.h),wt.nodes(nIdx{ii}));
    rangeOut = mat2cell(1:sum(rangeOut),1,rangeOut);
    [gtmp,atmp] = nodesMultid(nIdx{ii},rangeLoc,rangeOut,wt);
    g(end+1:end+numel(gtmp)) = gtmp;
    a(end+1:end+numel(atmp)) = atmp;
end
g = g(:);
a = a(:);


if nargout<2
   % Scale filters if a is not returned
   for nn=1:numel(g)
       g{nn}.h = g{nn}.h/sqrt(a(nn));
   end
end
















