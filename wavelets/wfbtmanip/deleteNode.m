function wt = deleteNode(nodeNo,wt)
%DELETENODE Removes specified node from the tree
%   Usage:  wt = deleteNode(nodeNo,wt)
%
%   Input parameters:
%         nodeNo   : Node index.
%         wt       : Structure containing description of the filter tree.
%
%   Output parameters:
%         wt       : Modified wt.

if(~isempty(find(wt.children{nodeNo}~=0)))
    error('Deleting non-leaf node!');
end

parId = wt.parents(nodeNo);
toZero = find(wt.children{parId}==nodeNo);
wt.children{parId}(toZero) = 0;

newIdx = 1:length(wt.nodes);
newIdx = newIdx(find(newIdx~=nodeNo));
wt.nodes = wt.nodes(newIdx);
%treeStruct.a = {treeStruct.a{newIdx}};
%treeStruct.origins = {treeStruct.origins{newIdx}};
wt.parents = wt.parents(newIdx); 
wt.children = wt.children(newIdx);

% and all children and parents with higher idx are lessened
 for ii =1:length(wt.children)
     biggerIdx = find(wt.children{ii}>nodeNo);
     wt.children{ii}(biggerIdx) = wt.children{ii}(biggerIdx)-1;
 end
 biggerIdx = find(wt.parents>nodeNo);
 wt.parents(biggerIdx) = wt.parents(biggerIdx)-1;