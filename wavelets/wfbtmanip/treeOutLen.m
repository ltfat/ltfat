function Lc = treeOutLen(L,doNoExt,wt)
%TREESUB

% All nodes with at least one final output.
termN = find(noOfNodeOutputs(1:numel(wt.nodes),wt)~=0);
% Range in filter outputs
outRange = rangeInLocalOutputs(termN,wt);
cRange = cell2mat(cellfun(@(rEl) rEl(:),rangeInOutputs(termN,wt),...
                  'UniformOutput',0));


Lctmp = nodeOutLen(termN,L,outRange,doNoExt,wt);
Lc = zeros(size(Lctmp));

Lc(cRange) = Lctmp;







