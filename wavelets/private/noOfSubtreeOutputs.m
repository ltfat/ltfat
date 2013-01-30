function noOut = noOfSubtreeOutputs(nodeNo,treeStruct)
noChildOut = noOfChildOutputs(nodeNo,treeStruct);
chan = max([length(treeStruct.nodes{nodeNo}.g), length(treeStruct.nodes{nodeNo}.h)]);
child = length(find(treeStruct.children{nodeNo}~=0));
noOut = chan -child + noChildOut;