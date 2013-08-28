function c=comp_wpfbt(f,wtNodes,rangeLoc,ext)
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
c = cell(sum(cellfun(@(wtEl) numel(wtEl.h),wtNodes)),1);

ca = f;
cOutRunIdx = 1;
cInRunIdxs = [1];
% Go over all nodes in breadth-first order
for jj=1:numel(wtNodes)
   % Node filters to a cell array
   hCell = cellfun(@(hEl) conj(flipud(hEl.h(:))),wtNodes{jj}.h(:),'UniformOutput',0);
   % Node filters subs. factors
   a = wtNodes{jj}.a;
   % Node filters initial skips
   if(doPer)
      offset = cellfun(@(hEl) 1-numel(hEl.h)-hEl.offset,wtNodes{jj}.h);
   else
      offset = -(a-1);
   end
   filtNo = numel(hCell);
   
   % Run filterbank
   c(cOutRunIdx:cOutRunIdx + filtNo-1)=...
                               comp_filterbank_td(ca,hCell,a,offset,ext);
   
   % Bookeeping. Store idxs of just computed outputs.
   outRange = cOutRunIdx:cOutRunIdx+filtNo-1;
   % Omit those, which are not decomposed further
   outRange(rangeLoc{jj}) = [];
   cInRunIdxs = [cInRunIdxs(2:end),outRange];
   
   cOutRunIdx = cOutRunIdx + filtNo;
   
   % Prepare input for the next iteration
   % Scaling introduced in order to preserve energy 
   % (parseval tight frame)
   % TO DO: Investigate how this scaling influences the |wpbest| 
   % algorithms.
   if ~isempty(cInRunIdxs)
      c{cInRunIdxs(1)} = c{cInRunIdxs(1)}/sqrt(2);
      ca = c{cInRunIdxs(1)};
   end
end   


