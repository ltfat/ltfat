function outRange = nodesLocOutRange(nodeNo,wt)
%NODESLOCOUTRANGE Node output index range of the terminal outputs
%   Usage:  outRange = nodesLocOutRange(nodeNo,wt);
%
%   Input parameters:
%         nodeNo     : Node index.
%         wt : Structure containing description of the filter tree.
%
%   Output parameters:
%         noOut      : Index range. 
%
%   `nodesLocOutRange(nodeNo,wt)` returns range of indexes of the
%   terminal outputs of the node `nodeNo`. For definition of the structure
%   see `wfbinit`.
%
%   See also: wfbtinit
%

if(any(nodeNo>numel(wt.nodes)))
   error('%s: Invalid node index range. Number of nodes is %d.\n',upper(mfilename),numel(wt.nodes));
end


nodesCount = length(nodeNo);
outRange = cell(nodesCount,1); 


nodeChans = cellfun(@(nEl) numel(nEl.g), wt.nodes(nodeNo));
chIdx = cellfun(@(chEl) find(chEl~=0), wt.children(nodeNo),'UniformOutput',0);


for ii = 1:nodesCount
 outRangeTmp = 1:nodeChans(ii);
 outRangeTmp(chIdx{ii}) = [];
 outRange{ii} = outRangeTmp;
end



