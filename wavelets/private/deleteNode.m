function treeStruct = deleteNode(nodeNo,treeStruct)
if(~isempty(find(treeStruct.children{nodeNo}~=0)))
    error('Deleting non-leaf node!');
end

parId = treeStruct.parents(nodeNo);
toZero = find(treeStruct.children{parId}==nodeNo);
treeStruct.children{parId}(toZero) = 0;

newIdx = 1:length(treeStruct.nodes);
newIdx = newIdx(find(newIdx~=nodeNo));
treeStruct.nodes = {treeStruct.nodes{newIdx}};
%treeStruct.a = {treeStruct.a{newIdx}};
%treeStruct.origins = {treeStruct.origins{newIdx}};
treeStruct.parents = treeStruct.parents(newIdx); 
treeStruct.children = {treeStruct.children{newIdx}}';

% and all children and parents with higher idx are lessened
 for ii =1:length(treeStruct.children)
     biggerIdx = find(treeStruct.children{ii}>nodeNo);
     treeStruct.children{ii}(biggerIdx) = treeStruct.children{ii}(biggerIdx)-1;
 end
 biggerIdx = find(treeStruct.parents>nodeNo);
 treeStruct.parents(biggerIdx) = treeStruct.parents(biggerIdx)-1;