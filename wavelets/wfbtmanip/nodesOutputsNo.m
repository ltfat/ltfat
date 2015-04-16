function noOut = nodesOutputsNo(nodeNo,wt)
%NODESOUTPUTSNO Number of node Outputs
%   Usage:  noOut = nodesOutputsNo(nodeNo,wt);
%
%   Input parameters:
%         nodeNo  : Node index.
%         wt      : Structure containing description of the filter tree.
%
%   Output parameters:
%         noOut      : Number of node outputs. 
%
%   `nodesOutputsNo(nodeNo,wt)` Return number of the terminal 
%   outputs of the node `nodeNo`. For definition of the structure
%   see `wfbinit`.
%
%   See also: wfbtinit
%

if(any(nodeNo>numel(wt.nodes)))
   error('%s: Invalid node index range. Number of nodes is %d.\n',upper(mfilename),numel(wt.nodes));
end

noOut = cellfun(@(nEl) numel(nEl.g), wt.nodes(nodeNo)) -...
        cellfun(@(chEl) numel(chEl(chEl~=0)), wt.children(nodeNo));




