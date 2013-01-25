function nodesIdxs = nodesBForder(treeStruct)
%NODESBFORDER Nodes in the Breadth-First search order
%  Usage:  nodesIdxs = nodesBForder(treeStruct)
%
%   Input parameters:
%         treeStruct  : Structure containing description of the filter tree.
%
%   Output parameters:
%         nodesIdxs   : Node indexes in the Breadth-First search order.
%
%   `nodesBForder(treeStruct)` For definition of the structure see
%   `wfbinit`.
%
%   See also: wfbtinit
%


%find root
nodeNo = find(treeStruct.parents==0);
nodesIdxs = [nodeNo];
toGoTrough = [nodeNo];

while ~isempty(toGoTrough)
   chtmp = find(treeStruct.children{toGoTrough(1)}~=0);
   chIdxtmp = treeStruct.children{toGoTrough(1)}(chtmp);
   nodesIdxs = [nodesIdxs,chIdxtmp];
   toGoTrough = [toGoTrough(2:end),chIdxtmp];
end