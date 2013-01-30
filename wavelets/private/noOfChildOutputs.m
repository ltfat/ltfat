function noOut = noOfChildOutputs(nodeNo,treeStruct)
noOut = 0;
childrenIdx = find(treeStruct.children{nodeNo}~=0);
children = treeStruct.children{nodeNo}(childrenIdx);
for nn=1:length(children)
   chNodeNo = children(nn);
   chan = max([length(treeStruct.nodes{chNodeNo}.g), length(treeStruct.nodes{chNodeNo}.h)]); 
   child = length(find(treeStruct.children{chNodeNo}~=0));
   noOut = noOut + chan -child;
   noOut = noOut + noOfChildOutputs(chNodeNo,treeStruct);
end