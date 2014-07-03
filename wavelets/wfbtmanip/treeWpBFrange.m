function [pOutIdxs,chOutIdxs] = treeWpBFrange(wt)

treePath = nodeBForder(0,wt);
trLen = numel(treePath);
pOutIdxs = zeros(1,trLen);
chOutIdxs = cell(1,trLen);
pRunIdx = [0];
chRunIdx = 1;
% do trough tree and look for nodeNo and its parent
for ii=1:trLen
    tmpfiltNo = length(wt.nodes{treePath(ii)}.g);
    locRange = nodesLocOutRange(treePath(ii),wt);
    diffRange = 1:tmpfiltNo;
    diffRange(locRange{1})=[];
    chOutIdxs{ii} = chRunIdx:chRunIdx+tmpfiltNo-1;
    chRunIdx = chRunIdx + tmpfiltNo;
    pOutIdxs(ii) = pRunIdx(1);
    pRunIdx = [pRunIdx(2:end),chOutIdxs{ii}(diffRange)];
end



pOutIdxs = pOutIdxs(end:-1:1); 
chOutIdxs = chOutIdxs(end:-1:1);


