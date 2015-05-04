function upsNo = nodesFiltUps(nodeNo,wt)
%NODEFILTUPS  Node upsamplig factor
%   Usage:  upsNo = nodesFiltUps(nodeNo,wt)
%
%   Input parameters:
%         wt  : Structure containing description of the filter tree.
%
%   Output parameters:
%         upsNo : Accumulated upsampling factor along path to root.
%
%   `nodesFiltUps(wt)` Returns upsampling factor, which can be used to
%   upsample the node filters using the a-trous algorithm.
%   For definition of the structure see |wfbinit|.
%
%   See also: wfbtinit
%

if(any(nodeNo>numel(wt.nodes)))
   error('%s: Invalid node index range. Number of nodes is %d.\n',upper(mfilename),numel(wt.nodes));
end

nodesCount = numel(nodeNo);
upsNo = zeros(nodesCount,1);
for ii=1:nodesCount
   tmpNodeNo = nodeNo(ii);
   upsNo(ii) = 1;
   while(wt.parents(tmpNodeNo))
       parentNo = wt.parents(tmpNodeNo);
       upsNo(ii)=upsNo(ii)*wt.nodes{parentNo}.a(wt.children{parentNo}==tmpNodeNo);
       tmpNodeNo=parentNo;
   end
end
