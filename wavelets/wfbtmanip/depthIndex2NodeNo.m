function [nodeNo,nodeChildIdx] = depthIndex2NodeNo(d,k,wt)
%DEPTHINDEX2NODENO Get node from depth and index in the tree
%   Usage: [nodeNo,nodeChildIdx] = depthIndex2NodeNo(d,k,wt)
%
%   `[nodeNo,nodeChildIdx] = depthIndex2NodeNo(d,k,wt)` returns node 
%   *nodeNo* and an array of its children nodes *nodeChildIdx* positioned
%   in depth *g* and index *k* in the tree *wt*.
%
if(d==0)
    nodeNo=0;
    nodeChildIdx=0;
    return;
end

% find ordered nodes at depth d-1
nodesNo = getNodesInDepth(d,wt);
if(isempty(nodesNo))
   error('%s: Depth of the tree is less than given d.',mfilename); 
end

% k is index in children of ordered nodes at depth d

nodeNo = zeros(numel(k),1);
nodeChildIdx = zeros(numel(k),1);
chNo = cumsum(cellfun( @(nEl) length(nEl.g),wt.nodes(nodesNo)));
chNoZ = [0;chNo(:)];

for kIdx=1:numel(k)
    ktmp = k(kIdx);
    idx = find(chNo>ktmp,1);
    if isempty(idx)
       error('%s: Index k=%i out of bounds.',mfilename,ktmp); 
    end    
    nodeNo(kIdx) = nodesNo(idx);
    nodeChildIdx(kIdx) = ktmp-chNoZ(idx)+1;
end

function nodd = getNodesInDepth(d,wt)
% find all nodes with d steps to the root ordered
if d==1
    % return root
    nodd = find(wt.parents==0);
    return;
end    

nbf = nodeBForder(0,wt);
nbfTmp = nbf;
tempd = 0;
while tempd<d
    nbf(nbfTmp==0) = [];
    nbfTmp(nbfTmp==0) = [];
    nbfTmp = wt.parents(nbfTmp);
    tempd = tempd+1;
end
nodd = nbf(nbfTmp==0);


