function demo_blockproc_basicloop(source,varargin)
%DEMO_BLOCKPROC_BASICLOOP Basic real-time audio manipulation
%


if nargin<1
   fprintf(['%s: To run the demo, use one of the following:\n',...
          'demo_blockproc_basicloop(''gspi.wav'') to play gspi.wav (any wav file will do).\n',...
          'demo_blockproc_basicloop(''dialog'') to choose the wav file via file chooser dialog GUI.\n',...
          'demo_blockproc_basicloop(f,''fs'',fs) to play from a column vector f using sampling frequency fs.\n',...
          'demo_blockproc_basicloop(''playrec'') to record from a mic and play simultaneously.\n',...
          'Avalable input and output devices can be listed by |blockdevices|.\n',...
          'Particular device can be chosen by passing additional key-value pair ''devid'',devid.\n',...
          'Output channels of the device cen be selected by additional key-value pair ''playch'',[ch1,ch2].\n',...
          'Input channels of the device cen be selected by additional key-value pair ''recch'',[ch1].\n',...
          ]...
          ,upper(mfilename));
    return;
end

try
   playrec('isInitialised');
catch
   error('%s: playrec or portaudio are not properly compiled. ',mfilename);
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