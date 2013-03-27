function f=comp_iuwpfbt(c,wtNodes,nodesUps,pOutIdxs,chOutIdxs)
%COMP_IUWPFBT Compute Inverse Undecimated Wavelet Packet Filter-Bank Tree
%   Usage:  f=comp_iuwpfbt(c,wtNodes,nodesUps,pOutIdxs,chOutIdxs)
%
%   Input parameters:
%         c          : Coefficients stored in L*M*W array.
%         wtNodes    : Filterbank tree nodes (elementary filterbans) in
%                      reverse BF order. Cell array of structures of length *nodeNo*.
%         nodesUps   : Filters upsampling factor of each node. Array of
%                      length *nodeNo*.
%         pOutIdxs   : Idx of each node's parent. Array of length *nodeNo*.
%         chOutIdxs  : Idxs of each node children. Cell array of vectors of
%                      length *nodeNo*.
%
%   Output parameters:
%         f     : Reconstructed data in L*W array.
%

% For each node in tree in the BF order...
 for jj=1:length(wtNodes)
    % Node filters subs. factors
    a = wtNodes{jj}.a;
    % Node filters to a matrix
    gMat = cell2mat(cellfun(@(gEl) gEl.h(:),wtNodes{jj}.filts(:)','UniformOutput',0));
    % Normalize each filter
    gMat = bsxfun(@rdivide,gMat,sqrt(a(:)'));
    % Node filters initial skips
    gDel = cellfun(@(gEl) gEl.d,wtNodes{jj}.filts);
    
    % Zero index position of the upsampled filters.
    skip = nodesUps(jj).*(gDel) - nodesUps(jj);

    % Run filterbank
    ctmp = comp_iatrousfilterbank_td(c(:,chOutIdxs{jj},:),gMat,nodesUps(jj),skip);
    
    if(pOutIdxs(jj))
       % Add to the existing subband
       c(:,pOutIdxs(jj),:) = 0.5*c(:,pOutIdxs(jj),:)+0.5*reshape(ctmp,size(ctmp,1),1,size(ctmp,2));
    else
       % We are at the root.
       f = ctmp;
    end
 end
