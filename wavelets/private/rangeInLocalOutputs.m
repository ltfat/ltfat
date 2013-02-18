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
nodesCount = length(nodeNo);
if(nodesCount>1)
   outRange = cell(nodesCount,1); 
else
   outRange = []; 
end

for ii = 1:nodesCount
   chIdx = find(treeStruct.children{nodeNo(ii)}~=0);
   chan = max([length(treeStruct.nodes{nodeNo(ii)}.g), length(treeStruct.nodes{nodeNo(ii)}.h)]);
   outNodes = zeros(chan,1);
   outNodes(chIdx) = treeStruct.children{nodeNo(ii)}(chIdx);
   if(iscell(outRange))
      outRange{ii} = find(outNodes==0);
   else
      outRange = find(outNodes==0);
   end
end

