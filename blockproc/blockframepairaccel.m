function [Fao,Fso] = blockframepairaccel(Fa, Fs, Lb, varargin)
%BLOCKACCEL Precompute structures for block processing
%   Usage: F = blockaccel(F,Lb);
%
%   `F=blockaccel(F,Lb)` have to be called for each frame object prior 
%   entering the main loop where |blockana| and |blocksyn| are called.
%   The function calls |frameaccel| and prepares structures for the
%   processing of a consecutive stream of blocks.
%
%
%
%      `'sliwin',sliwin`   : Slicing window. `sliwin` have to be a window
%                            of length *2L*. It is used in the slicing
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
      kv.anasliwin = 'rect';
   end

   Fao = blockframeaccel(Fa,Lb,'sliced','sliwin',kv.anasliwin);
   Fso = blockframeaccel(Fs,Lb,'sliced','sliwin',kv.synsliwin);
else
   Fao = blockframeaccel(Fa,Lb,flags.blockalg);
   Fso = blockframeaccel(Fs,Lb,flags.blockalg);
end




