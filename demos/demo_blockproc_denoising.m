function demo_blockproc_denoising(source,varargin)
%DEMO_BLOCKPROC_DENOISING Variable coefficients thresholding
%   Usage: demo_blockproc_denoising('gspi.wav')
%
%   For additional help call |demo_blockproc_denoising| without arguments.
%
%   The present demo allows you to set the coefficient threshold during the
%   playback using the control panel.
% 

if demo_blockproc_header(mfilename,nargin)
   return;
end

% Control pannel (Java object)
% Each entry determines one parameter to be changed during the main loop
% execution.
p = blockpanel({
               {'GdB','Gain',-20,20,0,21},...
               {'Thr','Treshold',0,0.1,0,1000}
               });
    
% Buffer length
bufLen = 1024;
% Number of frequency channels
M = 1000;

% Setup blocktream
fs=block(source,varargin{:},'loadind',p,'L',bufLen);

% Window length in ms
winLenms = 20; %floor(fs*winLenms/1e3)
[F,Fdual] = framepair('dgtreal',{'hann',floor(fs*winLenms/1e3)},'dual',40,M);
[Fa,Fs] = blockframepairaccel(F,Fdual, bufLen,'segola');


flag = 1;
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
   
  % Obtain parameters from the control panel
  gain = 10^(p.getParam('GdB')/20); % dB -> val
  thres = p.getParam('Thr');
  %bufLen = floor(p.getParam('bufLen'));

  % Read block of length bufLen
  [f,flag] = blockread();
  % Apply analysis frame
  c = blockana(Fa, f*gain); 
  % Plot
  % blockplot(fobj,F,c);
  % Apply thresholding
  c = thresh(c,thres,'soft');
  % Apply synthesis frame
  fhat = real(blocksyn(Fs, c, size(f,1)));
  % Play the block
  %fhat = f;
  blockplay(fhat);
end
blockdone(p);

