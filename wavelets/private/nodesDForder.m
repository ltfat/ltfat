function nodesIdxs = nodesDForder(treeStruct)
%NODESBFORDER Nodes in the Breadth-First search order
%  Usage:  nodesIdxs = nodesBForder(treeStruct)
%
%   Input parameters:
%         treeStruct  : Structure containing description of the filter tree.
%
%   Output parameters:
%         nodesIdxs   : Node indexes in the Breadth-First search order.
%
%   `nodesBForder(treeStruct)` For definition of the structure see
%   `wfbinit`.
%
%   See also: wfbtinit
%


%find root
nodeNo = find(treeStruct.parents==0);
nodesIdxs = [nodeNo,nodeSubtreeDF(nodeNo,treeStruct)];