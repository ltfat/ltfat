function subNo = nodeSub(nodeNo,wt)


if(any(nodeNo>numel(wt.nodes)))
   error('%s: Invalid node index range. Number of nodes is %d.\n',upper(mfilename),numel(wt.nodes));
end

nodeNoa = cellfun(@(nEl) nEl.a,wt.nodes,'UniformOutput',0);
nodeNoUps = nodeFiltUps(nodeNo,wt);

nodesCount = numel(nodeNo);
subNo = cell(1,nodesCount);
for ii=1:nodesCount
   subNo{ii} = nodeNoUps(ii)*nodeNoa{ii};
end

% if(nodesCount==1)
%    subNo = subNo{1};
% end


%subNo = nodeFiltUps(nodeNo,wt).*wt.nodes{nodeNo}.a;