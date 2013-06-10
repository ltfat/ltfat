function a = treeSub(wt)
%TREESUB

% All nodes with at least one final output.
termN = find(noOfNodeOutputs(1:numel(wt.nodes),wt)~=0);
% Range in filter outputs
outRangeTermN = rangeInLocalOutputs(termN,wt);

cRangeTermN = rangeInOutputs(termN,wt);

noOut = noOfOutputs(wt);
% Subsampling factors of the terminal nodes
subTermN = nodeSub(termN,wt);

a = zeros(noOut, 1);
for ii=1:numel(termN)
   a(cRangeTermN{ii}) = subTermN{ii}(outRangeTermN{ii});
end

