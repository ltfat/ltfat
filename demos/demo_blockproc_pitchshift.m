function demo_blockproc_pitchshift(source,varargin)
%DEMO_BLOCKPROC_PITCHSHIFT Pitch shift by Gabor coefficient bands shift
%   Usage: demo_blockproc_pitchshift('gspi.wav')
%
%   For additional help call |demo_blockproc_pitchshift| without arguments.
%
%   This script demonstrates a real-time Gabor coefficient manipulation.
%   Frequency bands are shifted up or down according to the slider
%   position.
%

if demo_blockproc_header(mfilename,nargin)
   return;
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
F = frame('dgtreal',{'hann',floor(fs*winLenms/1e3)},100,M);
Fdual = F;
[F,Fdual] = framepair('dgtreal',{'hann',floor(fs*winLenms/1e3)},'dual',128,M);
[Fa,Fs] = blockframepairaccel(F,Fdual, bufLen,'segola');

flag = 1;
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
   gain = blockpanelget(p,'GdB');
   gain = 10.^(gain/20);
   shift = fix(blockpanelget(p,'Shi'));

   [f,flag] = blockread();
   % Noise gate
   f=f*gain;
   
   c = blockana(Fa, f);
   
   % Actual coefficient shift
   cc = Fa.coef2native(c,size(c));
   
   
   if(strcmpi(source,'playrec'))
      % Hum removal (aka low-pass filter)
      cc(1:2,:,:) = 0;
   end
   
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
blockdone();
p.close();
fobj.close();
