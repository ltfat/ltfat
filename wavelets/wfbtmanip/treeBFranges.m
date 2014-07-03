function [nodesBF, rangeLoc, rangeOut] = treeBFranges(wt,varargin)

nodesBF = nodeBForder(0,wt);

if ~isempty(varargin(strcmp('rev',varargin)));
   nodesBF = fliplr(nodesBF); 
end

rangeLoc = nodesLocOutRange(nodesBF,wt);
rangeOut = nodesOutRange(nodesBF,wt);

