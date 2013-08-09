function demo_blockproc_pitchshift(source,varargin)
%DEMO_BLOCKPROC_PITCHSHIFT Pitch shift by Gabor coefficient bands shift
%
%   This script demonstrates a real-time Gabor coefficient manipulation.
%   Frequency bands are shifted up or down according to the slider
%   position.
%


if nargin<1
   fprintf(['\n%s: To run the demo, use one of the following:\n',...
          '%s(''gspi.wav'') to play gspi.wav (any wav file will do).\n',...
          '%s(''dialog'') to choose the wav file via file chooser dialog GUI.\n',...
          '%s(f,''fs'',fs) to play from a column vector f using sampling frequency fs.\n',...
          '%s(''playrec'') to record from a mic and play simultaneously.\n',...
          'Avalable input and output devices can be listed by |blockdevices|.\n',...
          'Particular device can be chosen by passing additional key-value pair ''devid'',devid.\n',...
          'Output channels of the device cen be selected by additional key-value pair ''playch'',[ch1,ch2].\n',...
          'Input channels of the device cen be selected by additional key-value pair ''recch'',[ch1].\n\n',...
          ]...
          ,upper(mfilename),mfilename,mfilename,mfilename,mfilename);
    return;
end

try
   playrec('isInitialised');
catch
   error('%s: playrec or portaudio are not properly compiled. ',mfilename);
end

M = 1000;

fobj = blockfigure();
            

% Basic Control pannel (Java object)
parg = {
        {'GdB','Gain',-20,20,0,21},...
        {'Shi','Shift',-200,200,0,401}
       };

p = blockpanel(parg);
            

bufLen = 1024;
% Setup blocktream
if isoctave
   fs=block(source,varargin{:},'L',bufLen);
else
   fs=block(source,varargin{:},'loadind',p,'L',bufLen);
end

% Window length in ms
winLenms = 20; 
[F,Fdual] = framepair('dgtreal',{'hann',floor(fs*winLenms/1e3)},'dual',100,M);
[Fa,Fs] = blockframepairaccel(F,Fdual, bufLen,'segola');

flag = 1;
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
   gain = blockpanelget(p,'GdB');
   gain = 10.^(gain/20);
   shift = fix(blockpanelget(p,'Shi'));

   [f,flag] = blockread();
   f=f*gain;
   
   c = blockana(Fa, f);
   
   % Actual coefficient shift
   cc = Fa.coef2native(c,size(c));
   if shift<0
      cc = [cc(-shift+1:end,:,:); zeros(-shift,size(cc,2),size(cc,3))];
   else
      cc = [zeros(shift,size(cc,2),size(cc,3)); cc(1:end-shift,:,:)];
   end
   c = Fa.native2coef(cc);
   
   blockplot(fobj,Fa,c(:,1));
   
   fhat = blocksyn(Fs, c, size(f,1));

   blockplay(fhat);
end
p.close();
fobj.close();
