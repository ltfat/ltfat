function nodesIdxs = nodeSubtreeBF(nodeNo,wt)
%NODESUBTREEBF Node subtree nodes in Breath-First order
%   Usage:  noOut = nodeSubtreeBF(nodeNo,wt);
%
%   Input parameters:
%         nodeNo  : Node index.
%         wt      : Structure containing description of the filter tree.
%
%   Output parameters:
%         noOut   : Nodes in a Breath-First order. 
%

% subtreeIdx = [];
% 
% children = treeStruct.children{nodeNo}(find(treeStruct.children{nodeNo}~=0));
% subtreeIdx(end+1:end+length(children)) = children;
% 
% for ii=1:length(children)
%    tmpSbIdx = nodeSubtreeBF(children(ii),treeStruct);
%    subtreeIdx(end+1:end+length(tmpSbIdx)) = tmpSbIdx;
% end

toGoTrough = [nodeNo];
nodesIdxs = [];
while ~isempty(toGoTrough)
   chtmp = find(wt.children{toGoTrough(1)}~=0);
   chIdxtmp = wt.children{toGoTrough(1)}(chtmp);
   nodesIdxs = [nodesIdxs,chIdxtmp];
   toGoTrough = [toGoTrough(2:end),chIdxtmp];
end