function c=comp_uwfbt(f,wtNodes,nodesUps,rangeLoc,rangeOut)
%COMP_WFBT Compute Undecimated Wavelet Filterbank Tree
%   Usage:  c=comp_uwfbt(f,wtNodes,nodesUps,rangeLoc,rangeOut);
%
%   Input parameters:
%         f        : Input L*W array.
%         wtNodes  : Filterbank tree nodes (elementary filterbanks) in
%                    BF order. Length *nodeNo* cell array of structures.
%         nodesUps : Filters upsampling factor of each node. Array of
%                    length *nodeNo*.
%         rangeLoc : Idxs of each node terminal outputs. Length *nodeNo* 
%                    cell array of vectors.
%         rangeOut : Output subband idxs of each node terminal outputs.
%
%   Output parameters:
%         c     : Coefficient array of dim. L*M*W.
%


% Pre-allocated output
[L, W] = size(f);
M = sum(cellfun(@(rEl) numel(rEl),rangeOut));
c = zeros(L,M,W);

% Convinience input reshape
ca = reshape(f,size(f,1),1,size(f,2));
% For each node in tree in the BF order...
for jj=1:numel(wtNodes)
   % Node filters subs. factors
   a = wtNodes{jj}.a;
   % Node filters to a matrix
   hMat = cell2mat(cellfun(@(hEl) hEl.h(:),wtNodes{jj}.filts(:)','UniformOutput',0));
   % Normalize each filter
   hMat = bsxfun(@rdivide,hMat,sqrt(a(:)'));
   % Node filters initial skips
   hDel = cellfun(@(hEl) hEl.d,wtNodes{jj}.filts);
   
   % Upsampling the filters.
   % hMatUps = comp_ups(hMat,nodesUps(jj),1);
   % Zero index position of the upsampled filters.
   skip = nodesUps(jj).*(hDel - 1);
   
   % Run filterbank.
   catmp=comp_atrousfilterbank_td(squeeze(ca(:,1,:)),hMat,nodesUps(jj),skip);
   % Bookkeeping
   % Copy what goes directly to the output...
   c(:,rangeOut{jj},:)=catmp(:,rangeLoc{jj},:);
   % ...and save the rest.
   diffRange = 1:size(hMat,2);
   diffRange(rangeLoc{jj}) = [];
   ca = horzcat(ca(:,2:end,:), catmp(:,diffRange,:));
end 


