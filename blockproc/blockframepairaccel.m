function [Fao,Fso] = blockframepairaccel(Fa, Fs, Lb, varargin)
%BLOCKFRAMEPAIRACCEL Precompute structures for block processing
%   Usage: F = blockframepairaccel(Fa,Fs,Lb);
%
%   `[Fao,Fso]=blockframepairaccel(Fa,Fs,Lb)` works similar to 
%   |blockframeaccel| with a pair of frames. The only difference from
%   calling |blockframeaccel| separatelly for each frame is correct
%   default choice of the slicing windows. Frame objects `Fa,Fs` will be
%   accelerated for length `2*Lb`.
%
%   The following optional arguments are recognized:
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
%      `'zpad',zpad`   : Number of zero samples the block will be padded
%                        after it is windowed by a slicing window. Note the
%                        frames will be accelerated for length
%                        `2*Lb+2*kv.zpad`. 
%

complainif_notenoughargs(nargin,3,'BLOCKFRAMEPAIRACCEL');
complainif_notvalidframeobj(Fa,'BLOCKFRAMEPAIRACCEL');
complainif_notvalidframeobj(Fs,'BLOCKFRAMEPAIRACCEL');

definput.flags.blockalg = {'naive','sliced','segola'};
definput.keyvals.anasliwin = [];
definput.keyvals.synsliwin = [];
definput.keyvals.zpad = 0;
[flags,kv]=ltfatarghelper({},definput,varargin);

isSliProp = ~isempty(kv.anasliwin) || ~isempty(kv.synsliwin) || kv.zpad~=0;
assert(~(~flags.do_sliced && isSliProp),...
   sprintf(['%s: Definig slicing window properties without setting the',...
            ' ''sliced'' flag.'], mfilename));

if flags.do_sliced 
   if isempty(kv.anasliwin)
      kv.anasliwin = 'hann';
   end
   
   if isempty(kv.synsliwin)
      kv.synsliwin = 'hann';
   end

   Fao = blockframeaccel(Fa,Lb,'sliced','sliwin',kv.anasliwin,...
                         'zpad',kv.zpad);
   Fso = blockframeaccel(Fs,Lb,'sliced','sliwin',kv.synsliwin,...
                         'zpad',kv.zpad);
else
   Fao = blockframeaccel(Fa,Lb,flags.blockalg);
   Fso = blockframeaccel(Fs,Lb,flags.blockalg);
end




