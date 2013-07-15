function wt = nat2freqOrder(wt)
%NAT2FREQORDER Natural To Frequency Ordering
%   Usage:  wt = nat2freqOrder(wt);
%
%   Input parameters:
%         wt    : Structure containing description of the filter tree.
%
%   Output parameters:
%         wt    : Structure containing description of the filter tree.
%
%   `nat2freqOrder(wt)` Creates new wavelet filterbank tree definition
%   with permuted order of some filters for purposes of the correct frequency
%   ordering of the resultant identical filters. For definition of the
%   structure see `wfbinit`.
%
%   See also: wfbtinit,  wfbtmultid, nodesBForder
%

treePath = nodesBForder(wt);
%skip root
treePath = treePath(2:end);

for ii=1:length(treePath)
    % should not be zero
    nodeId = treePath(ii);
    parentId = wt.parents(nodeId);
    % local index of the parent output connected to the treePath(ii) node,
    % is in range 1:chan
    locIdx = find(wt.children{parentId}==nodeId,1);
    
    % do nothing if the node is connected to the first (hopefully lowpass)
    % output
    if(rem(locIdx,2)~=1)
       % now for the filter reordering
       chan = numel(wt.nodes{nodeId}.g);
       wt.nodes{nodeId}.g = wt.nodes{nodeId}.g(chan:-1:1);
       wt.nodes{nodeId}.h = wt.nodes{nodeId}.h(chan:-1:1);
       wt.nodes{nodeId}.a = wt.nodes{nodeId}.a(chan:-1:1);
    end    
    
end


