function demo_blockproc_erblets(source,varargin)
%DEMO_BLOCKPROC_ERBLETS Basic real-time rolling erblet-spectrogram visualization
%   Usage: demo_blockproc_erblets('gspi.wav')
%
%   For additional help call |demo_blockproc_erblets| without arguments.
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
               });
            
fobj = blockfigure();
    
% Buffer length
% Larger the number the higher the processing delay. 1024 with fs=44100Hz
% makes ~23ms.
% The value can be any positive integer.
% Note that the processing itself can introduce additional delay.
bufLen = 1024;

% Setup blocktream
fs=block(source,varargin{:},'loadind',p,'L',bufLen);

% Number of filters
M = 200;
[g,a]=erbfilters(fs,'fractional','L',2*bufLen,'M',M,'real');
[F,Fdual] = framepair('filterbankreal','dual',g,a,size(a,1));
[Fa,Fs] = blockframepairaccel(F,Fdual,bufLen,'sliced');

flag = 1;
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
  % Get parameters 
  gain = blockpanelget(p,'GdB');
  gain = 10^(gain/20);

  % Read block of length bufLen
  [f,flag] = blockread();
  f = f*gain;
  % Apply analysis frame
  c = blockana(Fa, f); 
  % Plot
  blockplot(fobj,Fa,c(:,1));
  
  fhat = real(blocksyn(Fs, c, size(f,1)));

  blockplay(fhat);
end
blockdone(p,fobj);



