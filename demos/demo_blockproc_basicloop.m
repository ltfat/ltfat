function demo_blockproc_basicloop(source,varargin)
%DEMO_BLOCKPROC_BASICLOOP Basic real-time audio manipulation
%   Usage: demo_blockproc_basicloop('gspi.wav')
%
%   For additional help call |demo_blockproc_basicloop| without arguments.
%
%   The demo runs simple playback loop allowing to set gain in dB.
%

if demo_blockproc_header(mfilename,nargin)
   return;
end

% Basic Control pannel (Java object)
p = blockpanel({
               {'GdB','Gain',-20,20,0,21},...
               });

% Setup blocktream
block(source,varargin{:},'loadind',p);

flag = 1;
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
   gain = blockpanelget(p,'GdB');
   gain = 10^(gain/20);
   
   [f,flag] = blockread();
   blockplay(f*gain);
end
blockdone(p);
