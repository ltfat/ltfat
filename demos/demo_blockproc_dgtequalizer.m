function demo_blockproc_dgtequalizer(source,varargin)
%DEMO_BLOCKPROC_DGTEQUALIZER Basic real-time audio manipulation in the
%transform domain.
%
%   This script demonstrates a real-time Gabor coefficient manipulation.
%   Frequency bands of Gabor coefficients are multiplied (weighted) by
%   values taken from sliders having a similar effect as a octave equalizer.
%   The shown spectrogram is a result of a re-analysis of the synthetized 
%   block to show a frequency content of what is actually played. 
%


if nargin<1
   fprintf(['%s: To run the demo, use one of the following:\n',...
          'demo_blockproc_dgtequalizer(''gspi.wav'') to play gspi.wav (any wav file will do).\n',...
          'demo_blockproc_dgtequalizer(''dialog'') to choose the wav file via file chooser dialog GUI.\n',...
          'demo_blockproc_dgtequalizer(f,''fs'',fs) to play from a column vector f using sampling frequency fs.\n',...
          'demo_blockproc_dgtequalizer(''playrec'') to record from a mic and play simultaneously.\n',...
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

M = 1000;

fobj = blockfigure();
            
octaves = 6;
voices = 1;
eqbands = (octaves)*voices;

d = floor((floor(M/2)+1)*2.^(-(0:eqbands-1)./(voices)));
d = fliplr([d,0]); 

% Basic Control pannel (Java object)
parg = {{'GdB','Gain',-20,20,0,21}};
for ii=1:eqbands
   parg{end+1} = {sprintf('G%idB',ii),sprintf('band%i',ii),-20,20,0,21};
end

p = blockpanel(parg);
            

bufLen = 1024;
% Setup blocktream
if isoctave
   fs=block(source,varargin{:},'L',bufLen);
else
   fs=block(source,varargin{:},'loadind',p,'L',bufLen);
end

% Window length in ms
winLenms = 20; %floor(fs*winLenms/1e3)
[F,Fdual] = framepair('dgtreal',{'hann',floor(fs*winLenms/1e3)},'dual',200,M);
[Fa,Fs] = blockframepairaccel(F,Fdual, bufLen,'sliced');

flag = 1;
ola = [];
ola2 = [];
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
   gain = blockpanelget(p);
   gain = 10.^(gain/20);


   [f,flag] = blockread();
   f=f*gain(1);
   gain = gain(2:end);
   
   
   [c, ola] = blockana(Fa, f, ola);
   
   cc = framecoef2tf(Fa,c);
   % Do the weighting
   for ii=1:eqbands
      cc(d(ii)+1:d(ii+1),:,:) = gain(ii)*cc(d(ii)+1:d(ii+1),:,:);
   end
   c = frametf2coef(Fa,cc);
   
   fhat = blocksyn(Fs, c, size(f,1));
   
   
   blockplay(fhat);
   
   % Do re-analysis of the modified
   [c2, ola2] = blockana(Fa, fhat, ola2);
   blockplot(fobj,Fa,c2(:,1));
end
p.close();
fobj.close();
