function f = block_ifwt(c,w,J,Lb)
%BLOCK_IFWT IFWT wrapper for blockstream processing
%   Usage: f=block_ifwt(c,w,J,Lb);
%
%   `f = block_ifwt(c,w,J,Lb)` returns block of data reconstructed
%   from coefficients *c* using the SegDWT algorithm (based on overlap-add
%   block convolution) with wavelet filters *w* and *J* levels. The 
%   reconstructed block contains overlap to the next block(s).
%
%   Do not call this function directly. It is called from |blocksyn| when
%   using 'fwt' frame type with 'segola' block transform handling (see 
%   |blockframeaccel|).
%
%   Function should be independent of block_interface.
%
%   See also: block, block_ifwt
%
%   References: ltfatnote026

if nargin<4
   error('%s: Too few input parameters.',upper(mfilename));
end;


w = fwtinit(w);

%Lc = fwtclength(Lb,g,J,'per');
filtNo = length(w.g);
subbNo = (filtNo-1)*J+1;
Lc = zeros(subbNo,1);
runPtr = 0; 
levelLen = Lb;
  for jj=1:J
     for ff=filtNo:-1:2
        Lc(end-runPtr) = floor(levelLen/w.a(ff));
        runPtr = runPtr + 1;
     end
     levelLen = ceil(levelLen/w.a(1));
  end
Lc(1)=levelLen; 
c = mat2cell(c,Lc);

m = numel(w.g{1}.h);
a = w.a(1);
% Do the extension 
cstartZeros = zeros(numel(Lc),1);
filtNo = length(w.g);
runPtr = 0; 
for jj=1:J-1
   for ff=filtNo:-1:2
      cstartZeros(end-runPtr) = (a^(J-jj)-1)/(a-1)*(m-a);
      runPtr = runPtr + 1;
   end
end 

% Pad with zeros ()
cext = cellfun(@(cEl,cZEl) zeros(size(cEl,1)+cZEl,size(cEl,2)),c,num2cell(cstartZeros),'UniformOutput',0);
for jj=1:numel(cext)
   cext{jj}(end+1-size(c{jj},1):end,:) = c{jj};
end

Ls = Lb + (a^(J)-1)/(a-1)*(m-a);

%% ----- Run computation 
f = comp_ifwt(cext,w.g,w.a,J,Ls,'valid');
