function L = nodeInLen(nodeNo,inLen,doNoExt,treeStruct,varargin)
%NODEINLEN Length of the node input signal
%   Usage:  L = nodeInLen(nodeNo,inLen,doExt,treeStruct);
%
%   Input parameters:
%         nodeNo     : Node index.
%         inLen      : Filter thee input signal length.
%         doNoExt    : Expansive representation indicator.
%         treeStruct : Structure containing description of the filter tree.
%
%   Output parameters:
%         L          : Length of the node input signal 
%
%   `nodeInLen(nodeNo,inLen,doExt,treeStruct)` return length of the input
%   signal of the node `nodeNo`. For definition of the structure see `wfbinit`.
%
%   See also: wfbtinit
%

doSyn = 0;
if(~isempty(varargin))
    if(strcmpi(varargin{1},'syn'))
        doSyn = 1;
    end
end
subPat = [];
filtLenPat = [];
tmpNodeNo = nodeNo;

while(treeStruct.parents(tmpNodeNo))
   parentNo = treeStruct.parents(tmpNodeNo);
   tmpIdx = find(treeStruct.children{parentNo}==tmpNodeNo);
   subPat(end+1) = treeStruct.nodes{parentNo}.a(tmpIdx);
   filtLenPat(end+1) = length(treeStruct.nodes{parentNo}.filts{tmpIdx}.h);
%    if(doSyn)
%       filtLenPat(end+1) = length(treeStruct.nodes{parentNo}.g{tmpIdx});
%    else
%       filtLenPat(end+1) = length(treeStruct.nodes{parentNo}.h{tmpIdx}); 
%    end
   tmpNodeNo=parentNo;
end

subPat = subPat(end:-1:1);
filtLenPat = filtLenPat(end:-1:1);

L = inLen;
if(~doNoExt)
    for ii=1:length(subPat)
        L = floor((L+filtLenPat(ii)-1)/subPat(ii));
    end
else
    for ii=1:length(subPat)
        L = ceil(L/subPat(ii)); 
    end
end

