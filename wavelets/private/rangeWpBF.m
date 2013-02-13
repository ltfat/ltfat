function [pOutIdxs,chOutIdxs] = rangeWpBF(wt,varargin)

treePath = nodesBForder(wt);
trLen = length(treePath);
pOutIdxs = zeros(1,trLen);
chOutIdxs = cell(1,trLen);
pRunIdx = [0];
chRunIdx = 1;
% do trough tree and look for nodeNo and its parrent
for ii=1:trLen
    tmpfiltNo = length(wt.nodes{treePath(ii)}.filts);
    chOutIdxs{ii} = chRunIdx:chRunIdx+tmpfiltNo-1;
    chRunIdx = chRunIdx + tmpfiltNo;
    pOutIdxs(ii) = pRunIdx(1);
    pRunIdx = [pRunIdx(2:end),chOutIdxs{ii}];
end


if(~isempty(varargin))
    if(strcmpi(varargin{1},'rev'))
       pOutIdxs = pOutIdxs(end:-1:1); 
       chOutIdxs = chOutIdxs(end:-1:1);
    end
end

