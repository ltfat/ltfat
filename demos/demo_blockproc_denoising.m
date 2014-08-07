function demo_blockproc_denoising(source,varargin) %RUNASSCRIPT
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
    
        


% Number of frequency channels
M = 1000;

% Setup blocktream
try
    fs=block(source,varargin{:},'loadind',p);
catch
    % Close the windows if initialization fails
    blockdone(p);
    err = lasterror;
    error(err.message);
end

% Buffer length (30 ms)
bufLen = floor(30e-3*fs);

% Window length in ms
winLenms = 20; %floor(fs*winLenms/1e3)
[F,Fdual] = framepair('dgtreal',{'hann',floor(fs*winLenms/1e3)},'dual',40,M);
% Or using fwt
%[F,Fdual] = framepair('fwt','ana:symorth3','dual',7);
[Fa,Fs] = blockframepairaccel(F,Fdual, bufLen,'segola');


flag = 1;
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
   
  % Obtain parameters from the control panel
  gain = 10^(p.getParam('GdB')/20); % dB -> val
  thres = p.getParam('Thr');
  %bufLen = floor(p.getParam('bufLen'));

  % Read block of length bufLen
  [f,flag] = blockread(bufLen);
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

