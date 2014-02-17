function [Fao,Fso] = blockframepairaccel(Fa, Fs, Lb, varargin)
%BLOCKFRAMEPAIRACCEL Precompute structures for block processing
%   Usage: F = blockframepairaccel(Fa,Fs,Lb);
%
%   `[Fao,Fso]=blockframepairaccel(Fa,Fs,Lb)` works similar to 
%   |blockframeaccel| with a pair of frames. The only difference from
%   calling |blockframeaccel| separatelly for each frame is correct
%   default choice of the slicing windows.   
%
%      `'anasliwin',anasliwin`  : Analysis slicing window. `sliwin` have to
%                                 be a window of length *2Lb* or a string 
%                                 accepted by the |firwin| function. It is
%                                 used only in the slicing window approach.
%                                 The default is `'hann'`.
%
%      `'synsliwin',synsliwin`  : Synthesis slicing window. The same as the
%                                 previous one holds. The default is `'rect'`.
%
%      `'transa',transa`   : Transition area length. Number of zeros to be 
%                            used to pad the slicing window to *2Lb* from
%                            both sides. The option is used as padding
%                            only in the slicing window approach using 
%                            window defined as a string.
%


definput.flags.blockalg = {'naive','sliced','segola'};
definput.keyvals.anasliwin = [];
definput.keyvals.synsliwin = [];
definput.keyvals.transa = [];
[flags,kv]=ltfatarghelper({},definput,varargin);

assert(~(~flags.do_sliced && (~isempty(kv.anasliwin) || ...
         ~isempty(kv.synsliwin) || ~isempty(kv.transa))),...
   sprintf('%s: Definig slicing window without setting the ''silced'' flag.',...
   mfilename));

if flags.do_sliced 
   if isempty(kv.anasliwin)
      kv.anasliwin = 'hann';
   end
   
   if isempty(kv.synsliwin)
      kv.synsliwin = 'rect';
   end

   Fao = blockframeaccel(Fa,Lb,'sliced','sliwin',kv.anasliwin,...
                         'transa',kv.transa);
   Fso = blockframeaccel(Fs,Lb,'sliced','sliwin',kv.synsliwin,...
                         'transa',kv.transa);
else
   Fao = blockframeaccel(Fa,Lb,flags.blockalg);
   Fso = blockframeaccel(Fs,Lb,flags.blockalg);
end




