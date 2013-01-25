function outRange = rangeInLocalOutputs(nodeNo,treeStruct)
%RANGEINLOCALOUTPUTS Node output index range of the terminal outputs
%   Usage:  outRange = rangeInLocalOutputs(nodeNo,treeStruct);
%
%   Input parameters:
%         nodeNo     : Node index.
%         treeStruct : Structure containing description of the filter tree.
%
%   Output parameters:
%         noOut      : Index range. 
%
%   `rangeInLocalOutputs(nodeNo,treeStruct)` Return range of indexes of the
%   terminal outputs of the node `nodeNo`. For definition of the structure
%   see `wfbinit`.
%
%   See also: wfbtinit
%
chIdx = find(treeStruct.children{nodeNo}~=0);
chan = max([length(treeStruct.nodes{nodeNo}.g), length(treeStruct.nodes{nodeNo}.h)]);
outNodes = zeros(chan,1);
outNodes(chIdx) = treeStruct.children{nodeNo}(chIdx);
outRange = find(outNodes==0);


