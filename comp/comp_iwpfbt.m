function f=comp_iwpfbt(c,wtNodes,pOutIdxs,chOutIdxs,Ls,ext)
%COMP_IWFBT Compute Inverse Wavelet Packet Filter-Bank Tree
%   Usage:  f=comp_iwpfbt(c,wt,Ls,type,ext)
%
%   Input parameters:
%         c          : Coefficients stored in cell array.
%         wtNodes    : Filterbank tree nodes (elementary filterbans) in
%                      reverse BF order. Cell array of structures of length *nodeNo*.
%         pOutIdxs   : Idx of each node's parent. Array of length *nodeNo*.
%         chOutIdxs  : Idxs of each node children. Cell array of vectors of
%                      length *nodeNo*.
%         ext        : Type of the forward transform boundary handling.
%
%   Output parameters:
%         f          : Reconstructed data in L*W array.
%

% Do non-expansve transform if ext=='per'
doPer = strcmp(ext,'per');

% For each node in tree in the BF order...
 for jj=1:length(wtNodes)
    % Node filters to a cell array
    gCell = cellfun(@(gEl) gEl.h(:),wtNodes{jj}.filts(:),'UniformOutput',0);
    % Node filters subs. factors
    a = wtNodes{jj}.a;
    % Node filters initial skips
    if(doPer)
       skip = cellfun(@(gEl) gEl.d-1,wtNodes{jj}.filts);
    else
       skip = cellfun(@(gEl) numel(gEl),gCell) - 1 - (a - 1);
    end
    
    if(pOutIdxs(jj))
       % Run filterbank and add to the existing subband.
       ctmp = comp_ifilterbank_td(c(chOutIdxs{jj}),gCell,a,size(c{pOutIdxs(jj)},1),skip,ext);
       c{pOutIdxs(jj)} = 0.5*c{pOutIdxs(jj)}+0.5*ctmp;
    else
       % We are at the root.
       f = comp_ifilterbank_td(c(chOutIdxs{jj}),gCell,a,Ls,skip,ext);
    end
 end
     
 