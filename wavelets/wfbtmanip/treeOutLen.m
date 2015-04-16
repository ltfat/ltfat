function Lc = treeOutLen(L,doNoExt,wt)
%TREEOUTLEN  Lengths of tree subbands
%   Usage:  Lc = treeOutLen(L,doNoExt,wt)
%
%   Input parameters:
%         L       : Input signal length.
%         doNoExt : Flag. Expansive = false, Nonexpansive=true  
%         wt      : Structure containing description of the filter tree.
%
%   Output parameters:
%         Lc : Subband lengths.
%
%   `Lc = treeOutLen(L,doNoExt,wt)` returns lengths of tree subbands given
%   input signal length *L* and flag `doNoExt`. When true, the transform is
%   assumed to be non-expansive.
%   For definition of the structure see |wfbinit|.
%
%   See also: wfbtinit
%



[termN, rangeLoc, rangeOut] = treeBFranges(wt);
slice = ~cellfun(@isempty,rangeOut); % Limit to nodes with unconnected outputs
rangeOut=rangeOut(slice);
cRange = cell2mat(cellfun(@(rEl) rEl(:),rangeOut(:),...
                  'UniformOutput',0));

Lctmp = nodesOutLen(termN(slice),L,rangeLoc(slice),doNoExt,wt);
Lc = zeros(size(Lctmp));
Lc(cRange) = Lctmp;







