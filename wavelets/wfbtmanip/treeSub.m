function a = treeSub(wt)
%TREESUB  Identical subsampling factors
%   Usage:  a = treeSub(wt)
%
%   Input parameters:
%         wt  : Structure containing description of the filter tree.
%
%   Output parameters:
%         a : Subsampling factors.
%
%   `a = treeSub(wt)` returns subsampling factors asociated with the tree
%   subbands. For definition of the structure see |wfbinit|.
%
%   See also: wfbtinit
%

% Get nodes in BF order
nodesBF = nodeBForder(0,wt);
% All nodes with at least one final output.
termNslice = nodesOutputsNo(nodesBF,wt)~=0;
termN = nodesBF(termNslice);

% Range in filter outputs
outRangeTermN = nodesLocOutRange(termN,wt);

%cRangeTermN = nodesOutRange(termN,wt);
rangeOut = treeOutRange(wt);
% Get only nodes with some output
cRangeTermN = rangeOut(termNslice);

noOut = sum(cellfun(@numel,cRangeTermN));
% Subsampling factors of the terminal nodes
subTermN = nodesSub(termN,wt);

a = zeros(noOut, 1);
for ii=1:numel(termN)
   a(cRangeTermN{ii}) = subTermN{ii}(outRangeTermN{ii});
end

