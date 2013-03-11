function idxOut = idxOfNodeOutputs(nodeNo,wt)
%NOOFNODEOUTPUTS Number of node Outputs
%   Usage:  noOut = noOfNodeOutputs(nodeNo,treeStruct);
%
%   Input parameters:
%         nodeNo     : Node index.
%         treeStruct : Structure containing description of the filter tree.
%
%   Output parameters:
%         noOut      : Number of node outputs. 
%
%   `noOfNodeOutputs(nodeNo,treeStruct)` Return number of the terminal 
%   outputs of the node `nodeNo`. For definition of the structure
%   see `wfbinit`.
%
%   See also: wfbtinit
%

idxOut = find(wt.children{nodeNo}~=0);

