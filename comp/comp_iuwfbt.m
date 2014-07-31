function f=comp_iuwfbt(c,wtNodes,nodesUps,rangeLoc,rangeOut,scaling)
%COMP_IUWFBT Compute Inverse Undecimated Wavelet Filter-Bank Tree
%   Usage:  f=comp_iuwfbt(c,wtNodes,nodesUps,rangeLoc,rangeOut)
%
%   Input parameters:
%         c        : Coefficient array of dim. L*M*W.
%         wtNodes  : Filterbank tree nodes (elementary filterbanks) in
%                    BF order. Length *nodeNo* cell array of structures.
%         nodesUps : Filters upsampling factor of each node. Array of
%                    length *nodeNo*.
%         rangeLoc : Idxs of each node inputs. Length *nodeNo* 
%                    cell array of vectors.
%         rangeOut : Input subband idxs of each node inputs.
%
%   Output parameters:
%         f     : Reconstructed data L*W array.
%

L = size(c,1);
W = size(c,3);
catmp = [];
ca = [];
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
    offset = nodesUps(jj).*(gOffset) ;%- nodesUps(jj);
     
    % Re-allocate catmp if the filtNo differs from the one used in previous
    % iteration.
    filtNo = size(gMat,2);
    if(filtNo~=size(catmp,2))
       catmp = zeros(L,filtNo,W,class(c));
    end

    % Read from input subbands
    catmp(:,rangeLoc{jj},:) = c(:,rangeOut{jj},:);
    diffRange = 1:filtNo;
    diffRange(rangeLoc{jj}) = [];
    % Read from intermediate outputs
    if(~isempty(diffRange))
       catmp(:,diffRange(end:-1:1),:) = ca(:,1:numel(diffRange),:);
    end
    
    %Run filterbank
    catmp = comp_iatrousfilterbank_td(catmp,gMat,nodesUps(jj),offset);
    %Save intermediate output
    ca = horzcat(ca(:,numel(diffRange)+1:end,:),reshape(catmp,size(catmp,1),1,size(catmp,2)));
end
f = catmp;

