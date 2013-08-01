function demo_blockproc_slidingsgram(source,varargin)
%DEMO_BLOCKPROC_SLIDINGSGRAM Basic real-time rolling spectrogram visualization
%
% This demo shows a simple rolling spectrogram of whatever is specified in
% source. 


if nargin<1
   fprintf(['%s: To run the demo, use one of the following:\n',...
          'demo_blockproc_slidingsgram(''gspi.wav'') to play gspi.wav (any wav file will do).\n',...
          'demo_blockproc_slidingsgram(''dialog'') to choose the wav file via file chooser dialog GUI.\n',...
          'demo_blockproc_slidingsgram(f,''fs'',fs) to play from a column vector f using sampling frequency fs.\n',...
          'demo_blockproc_slidingsgram(''playrec'') to record from a mic and play simultaneously.\n',...
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

% Basic sepctrogram figure (Java object)
fobj = blockfigure();

%
bufLen = 1024;
% Setup blocktream
if isoctave
   fs=block(source,varargin{:},'L',bufLen);
else
   fs=block(source,varargin{:},'loadind',p,'L',bufLen);
end

% Using dgtreal with 20ms hann window, hop factor 80, 1000 channels.
% Redundancy factor 12.5
winLenms = 20;
a = 80;
M = 1000;

F = frame('dgtreal',{'hann',floor(fs*winLenms/1e3)},a,M);
F = blockframeaccel(F, bufLen,'segola');


flag = 1;
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
   % Obtain the global gain value
   gain = blockpanelget(p,'GdB');
   gain = 10^(gain/20);

   % Read the next block of samples
   [f,flag] = blockread();
   f=f*gain;
   
   % Do analysis using the specified frame. 
   c = blockana(F, f); 
   
   % Draw the first channel coefficients
   blockplot(fobj,F,c(:,1));
   
   % Enqueue to play
   blockplay(f);
end
p.close();
fobj.close();