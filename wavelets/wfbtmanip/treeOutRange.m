function outRange = treeOutRange(wt)
%TREEOUTRANGE Index range of the outputs
%   Usage:  outRange = treeOutRange(wt);
%
%   Input parameters:
%         wt         : Structure containing description of the filter tree.
%
%   Output parameters:
%         outRange   : Subband idx range.
%
%   `treeOutRange(nodeNo,wt)` returns index range in the global
%   tree subbands associated. For definition of the
%   structure see `wfbinit`.
%
%   See also: wfbtinit
%

nodesIdBF = nodeBForder(0,wt);
nodesIdBFinv = zeros(size(nodesIdBF));
nodesIdBFinv(nodesIdBF) = 1:numel(nodesIdBF);
% Filterbank structs in BF order
nodesFbBF = wt.nodes(nodesIdBF);
% Number of filters in each node in BF order
M = cellfun(@(nEl) numel(nEl.h),nodesFbBF);
% Number of unconnected outputs of nodes in BF order.
Munc = cellfun(@(nEl,chEl) numel(nEl.h)-numel(chEl(chEl~=0)),...
               nodesFbBF,wt.children(nodesIdBF));
           
childrenBForder = wt.children(nodesIdBF);

Munccumsum = [0,cumsum(Munc)];
outRange = zeros(1,sum(Munc));
idxLIFO = { {1,1,1} };
k = 1;

% n indexes are indices of BF order
while ~isempty(idxLIFO)
    % Pop first element
    [n,mstart,munc] = idxLIFO{end}{:}; idxLIFO = idxLIFO(1:end-1);
    for m=mstart:M(n)
        % If m is among children of the current node
        if m<=numel(childrenBForder{n}) && childrenBForder{n}(m) ~= 0 
            % Idex of next node in BF order
            nnext = nodesIdBFinv(childrenBForder{n}(m));
            if m<M(n)
               idxLIFO{end+1} = {n,m+1,munc}; 
            end
            idxLIFO{end+1} = {nnext,1,1};
            break;
        else
            outRange(munc + Munccumsum(n)) = k;
            munc = munc + 1;
            k = k+1;
        end
    end
end

outRange = mat2cell(outRange,1,Munc);





