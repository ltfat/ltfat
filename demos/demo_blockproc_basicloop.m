function demo_blockproc_basicloop(source,varargin) %RUNASSCRIPT
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
try
    fs = block(source,varargin{:},'loadind',p);
catch
    % Close the windows if initialization fails
    blockdone(p);
    err = lasterror;
    error(err.message);
end

% Set buffer length to 30 ms
L = floor(30e-3*fs);

flag = 1;
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
   gain = blockpanelget(p,'GdB');
   gain = 10^(gain/20);
   
   [f,flag] = blockread(L);
   % The following does nothing in the rec only mode.
   blockplay(f*gain);
   % The following does nothing if 'outfile' was not specified 
   blockwrite(f);
end
blockdone(p);
