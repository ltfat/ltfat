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

bufLen = 1024;
% Setup blocktream
if isoctave
   block(source,varargin{:},'L',bufLen);
else
   block(source,varargin{:},'loadind',p,'L',bufLen);
end

flag = 1;
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
   gain = blockpanelget(p,'GdB');
   gain = 10^(gain/20);
   
   [f,flag] = blockread();
   blockplay(f*gain);
end
p.close();