function [Fao,Fso] = blockframepairaccel(Fa, Fs, Lb, varargin)
%BLOCKFRAMEPAIRACCEL Precompute structures for block processing
%   Usage: F = blockframepairaccel(Fa,Fs,Lb);
%
%   `[Fao,Fso]=blockframepairaccel(Fa,Fs,Lb)` works similar to 
%   |blockframeaccel| with a pair of frames. The only difference from
%   calling |blockframeaccel| separatelly for each frame is correct
%   default choice of the slicing windows.   
%
%      `'sliwin',sliwin`   : Slicing window. `sliwin` have to be a window
%                            of length *2Lb*. It is used in the slicing
%                            window approach.

definput.flags.blockalg = {'naive','sliced','segola'};
definput.keyvals.anasliwin = [];
definput.keyvals.synsliwin = [];
[flags,kv]=ltfatarghelper({},definput,varargin);

assert(~(~flags.do_sliced && (~isempty(kv.anasliwin) || ~isempty(kv.synsliwin))),...
   '%s: Definig slicing window without setting the ''silced'' flag.',mfilename);

if flags.do_sliced 
   if isempty(kv.anasliwin)
      kv.anasliwin = 'hann';
   end
   
   if isempty(kv.synsliwin)
      kv.synsliwin = 'rect';
   end

   Fao = blockframeaccel(Fa,Lb,'sliced','sliwin',kv.anasliwin);
   Fso = blockframeaccel(Fs,Lb,'sliced','sliwin',kv.synsliwin);
else
   Fao = blockframeaccel(Fa,Lb,flags.blockalg);
   Fso = blockframeaccel(Fs,Lb,flags.blockalg);
end




