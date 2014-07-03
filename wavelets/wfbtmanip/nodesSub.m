function subNo = nodesSub(nodeNo,wt)


if(any(nodeNo>numel(wt.nodes)))
   error('%s: Invalid node index range. Number of nodes is %d.\n',upper(mfilename),numel(wt.nodes));
end

nodeNoa = cellfun(@(nEl) nEl.a,wt.nodes(nodeNo),'UniformOutput',0);
nodeNoUps = nodesFiltUps(nodeNo,wt);

nodesCount = numel(nodeNo);
subNo = cell(1,nodesCount);
for ii=1:nodesCount
   subNo{ii} = nodeNoUps(ii)*nodeNoa{ii};
end

