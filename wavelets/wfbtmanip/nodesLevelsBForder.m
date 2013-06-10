function nodesIdxs = nodesLevelsBForder(treeStruct)
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
toGoTrough = [nodeNo];
nodesIdxs = {nodeNo};
inLevel = [1];
counter = 0;
level = 2;
chIdxSum = 0;
while ~isempty(toGoTrough)
   chtmp = find(treeStruct.children{toGoTrough(1)}~=0);
   chIdxtmp = treeStruct.children{toGoTrough(1)}(chtmp);
   counter = counter + 1;

   if(length(nodesIdxs)<level&&~isempty(chIdxtmp))
       nodesIdxs = {nodesIdxs{:},[]}; 
   end
   
   chIdxSum = chIdxSum + length(chIdxtmp);
   if(~isempty(chIdxtmp))
       nodesIdxs{level} = [nodesIdxs{level},chIdxtmp];
   end
   
   toGoTrough = [toGoTrough(2:end),chIdxtmp];

   if(counter==inLevel(level-1))
       counter = 0;
       inLevel(level) = chIdxSum;
       level = level + 1;
       chIdxSum = 0;
   end
end

