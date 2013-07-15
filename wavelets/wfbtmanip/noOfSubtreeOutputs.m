function noOut = noOfSubtreeOutputs(nodeNo,wt)

noChildOut = noOfChildOutputs(nodeNo,wt);
chan = numel(wt.nodes{nodeNo}.g);
child = length(find(wt.children{nodeNo}~=0));
noOut = chan -child + noChildOut;