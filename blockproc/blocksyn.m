function fhat = blocksyn( F, c , Lb)
%BLOCKSYN Blockwise synthesis interface
%   Usage: blocksyn(F, c, Lb)
%
%   Input parameters:
%      F    : Frame object.
%      c    : Coefficients of a block.
%   Output parameters:
%      fhat : Reconstructed block.
%
%   `c=blocksyn(F, c, Lb)` reconstructs block `fhat` from coefficients `c`
%   using frame defined by `F`.

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(F)
  error('%s: First agument must be a frame definition structure.',upper(mfilename));
end;

% Next block index start (from a global point of view, starting with zero)
nextSb = block_interface('getPos');
% Block index start (from a global point of view, starting with zero)
Sb = nextSb-Lb;

switch(F.type)
   case 'fwt'
      J = F.J;
      w = F.g;
      m = numel(w.g{1}.h);
      a = w.a(1);
      blocksize = a^J;
      r = (a^J-1)/(a-1)*(m-1);
      Lbrec = (floor(nextSb/blocksize) - floor(Sb/blocksize))*blocksize;
      rSb = (a^J-1)/(a-1)*(m-a) + mod(Sb,a^J);
      over = r - rSb;
      f = block_ifwt(c,w,J,Lbrec);
      ol = loadOverlap(r-mod(Sb, a^J));
      olLen = size(ol,1);
      f(1:olLen-over,:) = f(1:olLen-over,:) + ol(1+over:end,:);
      f = [ol(1:over,:);f];
      storeOverlap(f,r-mod(nextSb, a^J));
      fhat = f(1:Lb,:);
   otherwise
      % General processing
      % Equal block length assumtion
      % Reconstruct
      f = frsyn(F,c);
      % Result should not be longer than 2*Lb
      f = f(1:2*Lb,:);
      % Load and add overlap (first half)
      ol = loadOverlap(Lb);
      olLen = size(ol,1);
      f(1:olLen,:) = f(1:olLen,:) + ol;
      % Store overlap (second half)
      storeOverlap(f,Lb);
      % Return first half
      fhat = f(1:Lb,:);
end

function overlap = loadOverlap(L)
   overlap = block_interface('getSynOverlap');
   if isempty(overlap)
     overlap = zeros(L,size(c,2));
   end
   Lo = size(overlap,1);
   if nargin<1
      L = Lo;
   end
   if L>Lo
      error('%s: Required more samples than stored.',upper(mfilename));
   end
   overlap = overlap(end-L+1:end,:);
end

function storeOverlap(fext,L)
   if L>size(fext)
       error('%s: Required more samples.',upper(mfilename));
   end
   block_interface('setSynOverlap',fext(end-L+1:end,:)); 
end
end