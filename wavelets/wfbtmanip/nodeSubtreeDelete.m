function wt = nodeSubtreeDelete(nodeNo,wt)
%DELETESUBTREE Removes subtree with root node
%   Usage:  wt = nodeSubtreeDelete(nodeNo,wt)
%
%   Input parameters:
%         nodeNo   : Node index.
%         wt       : Structure containing description of the filter tree.
%
%   Output parameters:
%         wt       : Modified wt.

complainif_notenoughargs(nargin,2,'DELETESUBTREE');
complainif_notposint(nodeNo,'DELETESUBTREE');

% All nodes to be deleted in breadth-first order
toDelete = nodeBForder(nodeNo,wt);

% Start deleting from the deepest nodes to avoid deleting nodes with
% children
for ii = length(toDelete):-1:1
  wt = nodeDelete(toDelete(ii),wt); 
  biggerIdx = toDelete>toDelete(ii);
  toDelete(biggerIdx) = toDelete(biggerIdx) - 1;
end

%wt = nodeDelete(nodeNo,wt); 

function wt = nodeDelete(nodeNo,wt)
%DELETENODE Removes specified node from the tree
%   Usage:  wt = nodeDelete(nodeNo,wt)
%
%   Input parameters:
%         nodeNo   : Node index.
%         wt       : Structure containing description of the filter tree.
%
%   Output parameters:
%         wt       : Modified wt.

complainif_notenoughargs(nargin,2,'DELETENODE');
complainif_notposint(nodeNo,'DELETENODE');

if any(wt.children{nodeNo}~=0)
    error('%s: Deleting a non-leaf node!',upper(mfilename));
end

% Removing a root node
if wt.parents(nodeNo)==0
    % Better clear all fields, than simply call wfbtinit
    fNames = fieldnames(wt);
    for ii=1:numel(fNames)
        wt.(fNames{ii})(:) = [];
    end
    return;
end

% Remove the node from it's parent children node list
parId = wt.parents(nodeNo);
wt.children{parId}(wt.children{parId}==nodeNo) = 0;

% newIdx = 1:length(wt.nodes);
% newIdx = newIdx(find(newIdx~=nodeNo));
% wt.nodes = wt.nodes(newIdx);
% wt.parents = wt.parents(newIdx); 
% wt.children = wt.children(newIdx);

% Remove the node from the structure completely
wt.nodes(nodeNo) = [];

if isfield(wt,'dualnodes')
    wt.dualnodes(nodeNo) = [];
end

wt.parents(nodeNo) = []; 
wt.children(nodeNo) = [];

% Since node was removed, the interconnections are now wrong. 
% Let's fix that.
for ii =1:length(wt.children)
    biggerIdx = wt.children{ii}>nodeNo;
    wt.children{ii}(biggerIdx) = wt.children{ii}(biggerIdx)-1;
end

% .. ant the same in the parents array
biggerIdx = wt.parents>nodeNo;
wt.parents(biggerIdx) = wt.parents(biggerIdx)-1;
