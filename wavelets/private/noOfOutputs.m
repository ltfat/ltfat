function noOut = noOfOutputs(wt)
%NOOFOUTPUTS Returns number of outputs of the filter tree 
%   Usage:  noOut=noOfOutputs(treeStruct)
%
%   Input parameters:
%         wt  : Structure containing description of the filter tree.
%
%   Output parameters:
%         noOut       : Number of outputs of the whole filter tree.
%
%   `noOfOutputs(wt)` For definition of the structure see `wfbinit`.
%
%   See also: wfbtinit
%

noOut = sum(noOfNodeOutputs(1:numel(wt.nodes),wt));

%noOut = sum( cellfun(@(nEl) numel(nEl.filts),wt.nodes) -...
%        cellfun(@(chEl) numel(chEl(chEl~=0)), wt.children) );

% Equivalent:     
% noOut = 0;
% for jj =1:length(wt.nodes)
%     chan = max([length(wt.nodes{jj}.filts), length(wt.nodes{jj}.h)]);
%     children = length(find(wt.children{jj}~=0));
%     noOut = noOut + chan-children;
% end