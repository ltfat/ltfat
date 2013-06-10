function wt = deleteSubtree(nodeNo,wt)
toDelete = nodeSubtreeBF(nodeNo,wt);

for ii = length(toDelete):-1:1
  wt = deleteNode(toDelete(ii),wt); 
  biggerIdx = find(toDelete>toDelete(ii));
  toDelete(biggerIdx) = toDelete(biggerIdx) - 1;
end
wt = deleteNode(nodeNo,wt); 