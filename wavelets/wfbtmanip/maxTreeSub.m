function a = maxTreeSub(wt)
%MAXTREESUB

% All nodes with at least one final output.
termN = find(noOfNodeOutputs(1:numel(wt.nodes),wt)~=0);
% Range in filter outputs
outRangeTermN = rangeInLocalOutputs(termN,wt);

% Subsampling factors of the terminal nodes
subTermN = nodeSub(termN,wt);

a = [];
for ii=1:numel(termN)
   a = [a; subTermN{ii}(outRangeTermN{ii})];
end

a = max(a);