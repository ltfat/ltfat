function demo_blockproc(source,varargin)
%DEMO_BLOCKPROC Real-time block processing demonstration
%   
%

if nargin<1
   error('%s: Too few input arguments.',upper(mfilename));
end

jo = javaObject('net.sourceforge.ltfat.ContFrame');
bufLen = 1024;
fs = block(source,'single',varargin{:});
%F = frame('dgtreal','gauss',10,100);
F = frame('fwt','sym10',4);
Fdual = framedual(F);
flag = 1;
while flag && jo.flag
  d = jo.shared;
  [f,flag] = blockread(bufLen);
  c = blockana(F, f); 
  c = thresh(c,d/100);
  fhat = blocksyn(Fdual, c, bufLen);
  blockplay(fhat);
end
jo.close();
