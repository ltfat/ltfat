function f=comp_iuwpfbt(c,wtNodes,nodesUps,pOutIdxs,chOutIdxs,scaling,interscaling)
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

interscalingfac = 1;
if strcmp('intscale',interscaling)
    interscalingfac = 1/2;
elseif strcmp('intsqrt',interscaling)
    interscalingfac = 1/sqrt(2);
end

% For each node in tree in the BF order...
 for jj=1:length(wtNodes)
    % Node filters subs. factors
    a = wtNodes{jj}.a;
    
    % Optionally scale the filters
    g = comp_filterbankscale(wtNodes{jj}.g(:),a(:),scaling);
    
    % Node filters to a matrix
    gMat = cell2mat(cellfun(@(gEl) gEl.h(:),g','UniformOutput',0));

    % Node filters initial skips
    gOffset = cellfun(@(gEl) gEl.offset,wtNodes{jj}.g);
    
    % Zero index position of the upsampled filters.
    offset = nodesUps(jj).*(gOffset);% + nodesUps(jj);

    % Run filterbank
    ctmp = comp_iatrousfilterbank_td(c(:,chOutIdxs{jj},:),gMat,nodesUps(jj),offset);
    
    if(pOutIdxs(jj))
       % Add to the existing subband
       c(:,pOutIdxs(jj),:) = interscalingfac*(c(:,pOutIdxs(jj),:)+reshape(ctmp,size(ctmp,1),1,size(ctmp,2)));
    else
       % We are at the root.
       f = ctmp;
    end
 end
