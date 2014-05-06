function demo_blockproc_slidingerblets(source,varargin) %RUNASSCRIPT
%DEMO_BLOCKPROC_SLIDINGERBLETS Basic real-time rolling erblet-spectrogram visualization
%   Usage: demo_blockproc_slidingerblets('gspi.wav')
%
%   For additional help call |demo_blockproc_slidingerblets| without arguments.
%
%   This demo shows a simple rolling erblet-spectrogram of whatever is specified in
%   source. 

if demo_blockproc_header(mfilename,nargin)
   return;
end

% Control pannel (Java object)
% Each entry determines one parameter to be changed during the main loop
% execution.
p = blockpanel({
               {'GdB','Gain',-20,20,0,21},...
               {'cMult','C mult',0,80,20,81}
               });
            
fobj = blockfigure();
    
% Buffer length
% Larger the number the higher the processing delay. 1024 with fs=44100Hz
% makes ~23ms.
% The value can be any positive integer.
% Note that the processing itself can introduce additional delay.
bufLen = 1024;
zpad = bufLen/2;

% Setup blocktream
fs=block(source,varargin{:},'loadind',p,'L',bufLen);

% Number of filters
M = 200;
F = frame('erbletfb',fs,2*bufLen+2*zpad,'fractionaluniform','M',M);
Fa = blockframeaccel(F,bufLen,'sliced','zpad',zpad);

flag = 1;
cola = [];
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
  % Get parameters 
  [gain,mult] = blockpanelget(p,'GdB','cMult');
  gain = 10^(gain/20);
  mult = 10^(mult/20);

  % Read block of length bufLen
  [f,flag] = blockread();
  f = f*gain;
  % Apply analysis frame
  c = blockana(Fa, f); 
  % Plot
  cola = blockplot(fobj,Fa,mult*c(:,1),cola);
  
  blockplay(f);
end
blockdone(p,fobj,Fa);



