function a = treeSub(wt)
%TREESUB

% All nodes with at least one final output.
termN = find(nodesOutputsNo(1:numel(wt.nodes),wt)~=0);
% Range in filter outputs
outRangeTermN = nodesLocOutRange(termN,wt);

cRangeTermN = nodesOutRange(termN,wt);

noOut = sum(cellfun(@numel,cRangeTermN));
% Subsampling factors of the terminal nodes
subTermN = nodesSub(termN,wt);

a = zeros(noOut, 1);
for ii=1:numel(termN)
   a(cRangeTermN{ii}) = subTermN{ii}(outRangeTermN{ii});
end

