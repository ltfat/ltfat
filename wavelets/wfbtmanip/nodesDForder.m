function nodesIdxs = nodesDForder(treeStruct,varargin)
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
%   Supported flags:
%
%             'ord','rev'
%
%   See also: wfbtinit
%


%find root
nodeNo = find(treeStruct.parents==0);
nodesIdxs = [nodeNo,nodeSubtreeDF(nodeNo,treeStruct)];


if(~isempty(varargin))
    if(strcmpi(varargin{1},'rev'))
       nodesIdxs = nodesIdxs(end:-1:1); 
    end
end