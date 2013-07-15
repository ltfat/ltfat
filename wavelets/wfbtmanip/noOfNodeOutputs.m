function noOut = noOfNodeOutputs(nodeNo,wt)
%NOOFNODEOUTPUTS Number of node Outputs
%   Usage:  noOut = noOfNodeOutputs(nodeNo,wt);
%
%   Input parameters:
%         nodeNo  : Node index.
%         wt      : Structure containing description of the filter tree.
%
%   Output parameters:
%         noOut      : Number of node outputs. 
%
%   `noOfNodeOutputs(nodeNo,wt)` Return number of the terminal 
%   outputs of the node `nodeNo`. For definition of the structure
%   see `wfbinit`.
%
%   See also: wfbtinit
%

if(any(nodeNo>numel(wt.nodes)))
   error('%s: Invalid node index range. Number of nodes is %d.\n',upper(mfilename),numel(wt.nodes));
end

%This is slow..
noOut = cellfun(@(nEl) numel(nEl.g), wt.nodes(nodeNo)) -...
        cellfun(@(chEl) numel(chEl(chEl~=0)), wt.children(nodeNo));

%This is even slower
% nodesCount = numel(nodeNo);
% noOut = zeros(nodesCount,1);
% for ii = 1:nodesCount
%    noOut(ii) = numel(wt.nodes{ii}.filts) - numel(find(wt.children{ii}~=0));
% end

%  chan = max([length(wt.nodes{nodeNo}.g), length(wt.nodes{nodeNo}.h)]);
%  child = length(find(wt.children{nodeNo}~=0));
%  noOut = chan -child;


