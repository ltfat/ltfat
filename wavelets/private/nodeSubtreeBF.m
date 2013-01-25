function subtreeIdx = nodeSubtreeBF(nodeNo,treeStruct)
subtreeIdx = [];

children = treeStruct.children{nodeNo}(find(treeStruct.children{nodeNo}~=0));
subtreeIdx(end+1:end+length(children)) = children;

for ii=1:length(children)
   tmpSbIdx = nodeSubtreeBF(children(ii),treeStruct);
   subtreeIdx(end+1:end+length(tmpSbIdx)) = tmpSbIdx;
end