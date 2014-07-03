function Lc = treeOutLen(L,doNoExt,wt)
%TREESUB

% All nodes with at least one final output.
termN = find(nodesOutputsNo(1:numel(wt.nodes),wt)~=0);
% Range in filter outputs
outRange = nodesLocOutRange(termN,wt);
cRange = cell2mat(cellfun(@(rEl) rEl(:),nodesOutRange(termN,wt),...
                  'UniformOutput',0));


Lctmp = nodesOutLen(termN,L,outRange,doNoExt,wt);
Lc = zeros(size(Lctmp));

Lc(cRange) = Lctmp;







