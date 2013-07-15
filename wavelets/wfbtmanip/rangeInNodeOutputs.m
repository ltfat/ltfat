function outRange = rangeInNodeOutputs(nodeNo,treeStruct)
%RANGEINNODEOUTPUTS Index range of the node outputs
%   Usage:  outRange = rangeInNodeOutputs(nodeNo,treeStruct)
%
%   Input parameters:
%         nodeNo     : Node index.
%         treeStruct : Structure containing description of the filter tree.
%
%   Output parameters:
%         outRange   : Index range. 
%
%   `rangeInNodeOutputs(nodeNo,treeStruct)` For definition of the structure
%   see `wfbinit`.
%
%   See also: wfbtinit
%
chIdx = find(treeStruct.children{nodeNo}~=0);
chan = numel(treeStruct.nodes{nodeNo}.g);
outNodes = zeros(chan,1);
outNodes(chIdx) = treeStruct.children{nodeNo}(chIdx);

outRangeStart = 0;
outRange = [];
for ii=1:chan
   if(outNodes(ii)==0)
       outRange(end+1) = outRangeStart+1;
       outRangeStart = outRangeStart+1;
   else
      outRangeStart=outRangeStart + noOfSubtreeOutputs(outNodes(ii),treeStruct);
   end
end