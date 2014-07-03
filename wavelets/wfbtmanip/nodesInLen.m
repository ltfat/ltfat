function L = nodesInLen(nodeNo,inLen,doNoExt,wt)
%NODESINLEN Length of the node input signal
%   Usage:  L = nodesInLen(nodeNo,inLen,doExt,treeStruct);
%
%   Input parameters:
%         nodeNo     : Node index.
%         inLen      : Filter thee input signal length.
%         doNoExt    : Expansive representation indicator.
%         wt         : Structure containing description of the filter tree.
%
%   Output parameters:
%         Lin        : Length of the node input signal 
%
%   `nodesInLen(nodeNo,inLen,doExt,treeStruct)` return length of the input
%   signal of the node `nodeNo`. For definition of the structure see `wfbinit`.
%
%   See also: wfbtinit
%

L = zeros(numel(nodeNo),1);
for nn=1:length(nodeNo)
    subPat = [];
    filtLenPat = [];
    tmpNodeNo = nodeNo(nn);

    while(wt.parents(tmpNodeNo))
       parentNo = wt.parents(tmpNodeNo);
       tmpIdx = find(wt.children{parentNo}==tmpNodeNo);
       subPat(end+1) = wt.nodes{parentNo}.a(tmpIdx);
       filtLenPat(end+1) = length(wt.nodes{parentNo}.g{tmpIdx}.h);
       tmpNodeNo=parentNo;
    end

    subPat = subPat(end:-1:1);
    filtLenPat = filtLenPat(end:-1:1);

    L(nn) = inLen;
    if(~doNoExt)
        for ii=1:length(subPat)
            L(nn) = floor((L(nn)+filtLenPat(ii)-1)/subPat(ii));
        end
    else
        for ii=1:length(subPat)
            L(nn) = ceil(L(nn)/subPat(ii)); 
        end
    end
end

