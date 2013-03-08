function nodesIdxs = nodeSubtreeDF(nodeNo,wt)
%NODESUBTREEBF Node subtree nodes in Depth-First order
%   Usage:  noOut = nodeSubtreeDF(nodeNo,wt);
%
%   Input parameters:
%         nodeNo  : Node index.
%         wt      : Structure containing description of the filter tree.
%
%   Output parameters:
%         noOut   : Nodes in a Depth-First order. 
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

toGoTrough = nodeNo;
nodesIdxs = [];
while ~isempty(toGoTrough)
   chtmp = find(wt.children{toGoTrough(1)}~=0);
   chIdxtmp = wt.children{toGoTrough(1)}(chtmp);
   nodesIdxs = [nodesIdxs,toGoTrough(1)];
   toGoTrough = [chIdxtmp,toGoTrough(2:end)];
end

% remove the nodeNo. Just to be consistent with nodeSubtreeBF
% TO DO: is it wise?
nodesIdxs = nodesIdxs(2:end);