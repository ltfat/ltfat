function noOut = noOfSubtreeOutputs(nodeNo,wt)

noChildOut = noOfChildOutputs(nodeNo,wt);
chan = max([length(wt.nodes{nodeNo}.g), length(wt.nodes{nodeNo}.h)]);
child = length(find(wt.children{nodeNo}~=0));
noOut = chan -child + noChildOut;