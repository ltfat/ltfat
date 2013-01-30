function treeStruct = deleteSubtree(nodeNo,treeStruct)
toDelete = nodeSubtreeBF(nodeNo,treeStruct);

for ii = length(toDelete):-1:1
  treeStruct = deleteNode(toDelete(ii),treeStruct); 
  biggerIdx = find(toDelete>toDelete(ii));
  toDelete(biggerIdx) = toDelete(biggerIdx) - 1;
end
treeStruct = deleteNode(nodeNo,treeStruct); 