function Lc = nodesOutLen(nodeNo,L,outRange,doNoExt,wt)
%NODESOUTLEN Length of the node output
%   Usage:  Lc = nodesOutLen(nodeNo,inLen,doExt,wt);
%
%   Input parameters:
%         nodeNo     : Node index(es).
%         inLen      : Filter thee input signal length.
%         outRange   : Cell array. Each element is a vector of local out.
%         indexes.
%         doNoExt    : Expansive representation indicator.
%         wt         : Structure containing description of the filter tree.
%
%   Output parameters:
%         Lin        : Length of the node input signal 
%
%   `nodesOutLen(nodeNo,inLen,doExt,treeStruct)` return length of the input
%   signal of the node `nodeNo`. For definition of the structure see `wfbinit`.
%
%   See also: wfbtinit
%
if isempty(outRange)
    outRange = cellfun(@(nEl) 1:numel(nEl.g),wt.nodes(nodeNo),'UniformOutput',0);
end

Lc = zeros(sum(cellfun(@numel,outRange)),1);

inLens = nodesInLen(nodeNo,L,doNoExt,wt);

Lcidx = 1;
for ii=1:numel(inLens)
    nodeHlen = cellfun(@(nEl) numel(nEl.h),...
               wt.nodes{nodeNo(ii)}.g(outRange{ii}));
    nodea =  wt.nodes{nodeNo(ii)}.a(outRange{ii});
    
    if(~doNoExt)
       Lc(Lcidx:Lcidx+numel(nodeHlen)-1) = floor((inLens(ii)...
                                            +nodeHlen(:)-1)./nodea(:));
    else
       Lc(Lcidx:Lcidx+numel(nodeHlen)-1) = ceil(inLens(ii)./nodea(:)); 
    end
    Lcidx = Lcidx + numel(nodeHlen);
end





