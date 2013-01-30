function noOut = noOfOutputs(treeStruct)
%NOOFOUTPUTS Returns number of outputs of the filter tree 
%   Usage:  noOut=noOfOutputs(treeStruct)
%
%   Input parameters:
%         treeStruct  : Structure containing description of the filter tree.
%
%   Output parameters:
%         noOut       : Number of outputs of the filter tree.
%
%   `noOfOutputs(treeStruct)` For definition of the structure see
%   `wfbinit`.
%
%   See also: wfbtinit
%
noOut = 0;
for jj =1:length(treeStruct.nodes)
    chan = max([length(treeStruct.nodes{jj}.g), length(treeStruct.nodes{jj}.h)]);
    children = length(find(treeStruct.children{jj}~=0));
    noOut = noOut + chan-children;
end