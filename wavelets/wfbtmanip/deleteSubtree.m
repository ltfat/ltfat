function wt = deleteSubtree(nodeNo,wt)
%DELETESUBTREE Removes subtree with root node
%   Usage:  wt = deleteSubtree(nodeNo,wt)
%
%   Input parameters:
%         nodeNo   : Node index.
%         wt       : Structure containing description of the filter tree.
%
%   Output parameters:
%         wt       : Modified wt.

toDelete = nodeSubtreeBF(nodeNo,wt);

for ii = length(toDelete):-1:1
  wt = deleteNode(toDelete(ii),wt); 
  biggerIdx = find(toDelete>toDelete(ii));
  toDelete(biggerIdx) = toDelete(biggerIdx) - 1;
end
wt = deleteNode(nodeNo,wt); 