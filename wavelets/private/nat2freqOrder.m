function treeStruct = nat2freqOrder(treeStruct)
%NAT2FREQORDER Natural To Frequency Ordering
%   Usage:  wtree = nat2freqOrder(doSyn,wtree);
%
%   Input parameters:
%         doSyn    : 0 - analysis part, 1 - synthesis part, 2 or more
%         both
%         wtree    : Structure containing description of the filter tree.
%
%   Output parameters:
%         wtree    : Structure containing description of the filter tree.
%
%   `nat2freqOrder(doSyn,wtree)` Creates new wavelet filterbank tree definition
%   with permuted order of some filters for purposes of the correct frequency
%   ordering of the resultant identical filters. For definition of the
%   structure see `wfbinit`.
%
%   See also: wfbtinit,  wfbtmultid, nodesBForder
%

treePath = nodesBForder(treeStruct);
%skip root
treePath = treePath(2:end);

for ii=1:length(treePath)
    % should not be zero
    nodeId = treePath(ii);
    parentId = treeStruct.parents(nodeId);
    % number of parent outputs
    chan = length(treeStruct.nodes{parentId}.filts);
%     if(doSyn)
%        chan = length(treeStruct.nodes{parentId}.g);
%     else
%        chan = length(treeStruct.nodes{parentId}.h);
%     end
    % local index of the parent output connected to the treePath(ii) node,
    % is in range 1:chan
    locIdx = find(treeStruct.children{parentId}==nodeId,1);
    
    % do nothing if the node is connected to the first (hopefully lowpass)
    % output
    if(rem(locIdx,2)~=1)
       % now for the filter reordering
       range = chan:-1:1;
       treeStruct.nodes{nodeId}.filts = {treeStruct.nodes{nodeId}.filts{range}};
%        if(doSyn)
%           treeStruct.nodes{nodeId}.g = {treeStruct.nodes{nodeId}.g{range}};
%        else
%           treeStruct.nodes{nodeId}.h = {treeStruct.nodes{nodeId}.h{range}};
%        end
    end    
    
end


