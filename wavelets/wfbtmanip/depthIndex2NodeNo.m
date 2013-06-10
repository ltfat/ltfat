function [nodeNo,nodeChildIdx] = depthIndex2NodeNo(d,k,treeStruct)

if(d==0)
    nodeNo=0;
    nodeChildIdx=0;
    return;
end

% find ordered nodes at depth d-1
nodesNo = getNodesInDepth(d,treeStruct);
if(isempty(nodesNo))
   error('%s: Depth of the tree is less than given d.',mfilename); 
end
ktemp = k;
% k is index in children of ordered nodes at depth d
for ii=1:length(nodesNo)
    chNo = max([length(treeStruct.nodes{nodesNo(ii)}.g), length(treeStruct.nodes{nodesNo(ii)}.h)]);
    %chNo = length(treeStruct.nodes{nodesNo(ii)});
    if(ktemp<chNo)
        nodeChildIdx = ktemp+1;
        nodeNo = nodesNo(ii);
        return;
    else
        ktemp = ktemp-chNo;
    end
end

error('%s: Index k out of bounds.',mfilename);


function nodd = getNodesInDepth(d,treeStruct)
% find all nodes with d steps to the root ordered
nodd = [];
toGoTrough = {};


   nodeNo = find(treeStruct.parents==0);
   toGoTrough = cell(d+1,1);
   toGoTrough{1} = nodeNo;
   tempd = 1;


while(tempd<d)
    for jj=1:length(toGoTrough{tempd})
       actNod = toGoTrough{tempd}(jj);
       childrenIdx = find(treeStruct.children{actNod}~=0);
       ch = treeStruct.children{actNod}(childrenIdx);
       toGoTrough{tempd+1} = [toGoTrough{tempd+1},ch];
    end

    tempd = tempd+1;
end

nodd=toGoTrough{d};
