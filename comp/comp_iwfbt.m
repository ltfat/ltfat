function f=comp_iwfbt(c,wtNodes,outLens,rangeLoc,rangeOut,ext)
%COMP_IWFBT Compute Inverse Wavelet Filter-Bank Tree
%   Usage:  f=comp_iwfbt(c,wtNodes,outLens,rangeLoc,rangeOut,ext)
%
%   Input parameters:
%         c        : Coefficients stored in the cell array.
%         wtNodes  : Filterbank tree nodes (elementary filterbanks) in
%                    reverse BF order. Length *nodeNo* cell array of structures.
%         outLens  : Output lengths of each node. Length *nodeNo* array.
%         rangeLoc : Idxs of each node inputs. Length *nodeNo* 
%                    cell array of vectors.
%         rangeOut : Input subband idxs of each node inputs.
%         ext      : Type of the forward transform boundary handling.
%
%   Output parameters:
%         f       : Reconstructed outLens(end)*W array.
%

% Do non-expansve transform if ext='per'
doPer = strcmp(ext,'per');

ca = {};
 % Go over all nodes in breadth-first order
for jj=1:length(wtNodes)
    % Node filters to a cell array
    % gCell = cellfun(@(gEl) conj(flipud(gEl.h(:))),wtNodes{jj}.g(:),'UniformOutput',0);
    gCell = cellfun(@(gEl) gEl.h(:),wtNodes{jj}.g(:),'UniformOutput',0);
    % Node filters subs. factors
    a = wtNodes{jj}.a;
    % Node filters initial skips
    if(doPer)
       % offset = cellfun(@(gEl) 1-numel(gEl.h)-gEl.offset,wtNodes{jj}.g(:));
       offset = cellfun(@(gEl) gEl.offset,wtNodes{jj}.g(:));
    else
       offset = -(a-1);
    end
    filtNo = numel(gCell);
    
    % Prepare input cell-array
    catmp = cell(filtNo,1);
    % Read data from subbands
    catmp(rangeLoc{jj}) = c(rangeOut{jj});
    diffRange = 1:filtNo;
    diffRange(rangeLoc{jj}) = [];
    % Read data from intermediate outputs (filters are taken in reverse order)
    catmp(diffRange(end:-1:1)) = ca(1:numel(diffRange));
    
    %Run filterbank
    catmp = comp_ifilterbank_td(catmp,gCell,a,outLens(jj),offset,ext);
    %Save intermediate output
    ca = [ca(numel(diffRange)+1:end);catmp];
end
f = catmp;

