function noOut = noOfChildOutputs(nodeNo,wt)


noOut = 0;
childrenIdx = find(wt.children{nodeNo}~=0);
children = wt.children{nodeNo}(childrenIdx);
for nn=1:length(children)
   chNodeNo = children(nn);
   chan = numel(wt.nodes{chNodeNo}.g); 
   child = numel(find(wt.children{chNodeNo}~=0));
   noOut = noOut + chan - child;
   noOut = noOut + noOfChildOutputs(chNodeNo,wt);
end