function c=comp_wpfbt(f,wtNodes,ext)
%COMP_WPFBT Compute Wavelet Packet Filterbank Tree
%   Usage:  c=comp_wpfbt(f,wtNodes,ext);
%
%   Input parameters:
%         f        : Input L*W array.
%         wtNodes  : Filterbank tree nodes (elementary filterbanks) in
%                    BF order. Length *nodeNo* cell array of structures.
%         ext      : Type of the forward transform boundary handling.
%
%   Output parameters:
%         c        : Coefficients stored in cell-array. Each element is one
%                    subband (matrix with W columns).
%

% Do non-expansve transform if ext=='per'
doPer = strcmp(ext,'per');
% Pre-allocated output
c = cell(sum(cellfun(@(wtEl) numel(wtEl.filts),wtNodes)),1);

ca = f;
cOutRunIdx = 1;
cInRunIdxs = [1];
% Go over all nodes in breadth-first order
for jj=1:numel(wtNodes)
   % Node filters to a cell array
   hCell = cellfun(@(hEl) hEl.h(:),wtNodes{jj}.filts(:),'UniformOutput',0);
   % Node filters subs. factors
   a = wtNodes{jj}.a;
   % Node filters initial skips
   if(doPer)
      skip = cellfun(@(hEl) hEl.d-1,wtNodes{jj}.filts);
   else
      skip = a-1;
   end
   filtNo = numel(hCell);

   % Run filterbank
   c(cOutRunIdx:cOutRunIdx + filtNo-1)=...
                               comp_filterbank_td(ca,hCell,a,skip,ext);
   
   % Bookeeping. Store idxs of just computed outputs.
   cInRunIdxs = [cInRunIdxs(2:end),cOutRunIdx:cOutRunIdx+filtNo-1];
   cOutRunIdx = cOutRunIdx + filtNo;
   
   % Prepare input for the next iteration
   ca = c{cInRunIdxs(1)};
end   


