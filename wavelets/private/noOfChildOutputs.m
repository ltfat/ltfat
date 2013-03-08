function noOut = noOfChildOutputs(nodeNo,wt)
noOut = 0;
childrenIdx = find(wt.children{nodeNo}~=0);
children = wt.children{nodeNo}(childrenIdx);
for nn=1:length(children)
   chNodeNo = children(nn);
   chan = max([length(wt.nodes{chNodeNo}.g), length(wt.nodes{chNodeNo}.h)]); 
   child = length(find(wt.children{chNodeNo}~=0));
   noOut = noOut + chan -child;
   noOut = noOut + noOfChildOutputs(chNodeNo,wt);
end