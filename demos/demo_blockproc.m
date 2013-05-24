function demo_blockproc(source,varargin)
%DEMO_BLOCKPROC Real-time block processing demonstration
% Usage: demo_blockproc('gspi.wav')
%        demo_blockproc('playrec')
%   
% The demo shows capabilities of the LTFAT block-stream processing
% framework which makes it possible (to some extend) to do the processing
% in a real time during the audio playback or recording.
%
% The main loop consists of several steps: It reads block of samples from
% the input (wav file or mic) and analyses it using the defined frame. The
% coefficients are thresholded before they are used for the synthesis using
% the dual frame. Finally, the processed block is sent to the sound output.
% With a suitable frame, such processing works particularly well for a
% background noise reduction. To fully evaluate the denoising, you can
% increase the background noise artificially by increasing the sensitivity
% of your microphone.
%
% The present demo allows you to set the coefficient threshold during the
% playback using the control panel.
%
% PLEASE NOTE: Matlab is far from being real-time processing friendly so
% the playback can skip when you make your system to do something else. The
% playback quality could be improved by setting the Matlab process a higher
% run priority.

if nargin<1
   error(['%s: To run the demo, use one of the following:\n',...
          'demo_blockproc(''gspi.wav'') to play gspi.wav (any wav file will do).\n',...
          'demo_blockproc(''playrec'') to record from a mic and play processed simultaneously.']...
          ,upper(mfilename));
end

% Control pannel (Java object)
% Each entry determines one parameter to be changed during the main loop
% execution.
p = blockpanel({
               {'GdB','Gain',-20,20,0,21},...
               {'Thr','Treshold',0,0.1,0,1000},
               });
    
% Buffer length
% Larger the number the higher the processing delay. 1024 with fs=44100Hz
% makes ~23ms.
% The value can be any positive integer.
% Note that the processing itself can introduce additional delay.
bufLen = 1024;

% Setup blocktream
fs = block(source,'nbuf',1,'single',varargin{:});

% Choose a frame and contruct the dual
%F = frame('dgtreal','gauss',10,100);
F = frame('fwt','db8',1);
Fdual = framedual(F);

flag = 1;
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag
  % Obtain parameters from the control panel
  gain = 10^(p.getParam('GdB')/20); % dB -> val
  thres = p.getParam('Thr');
  % Read block of length bufLen
  [f,flag] = blockread(bufLen);
  % Apply gain
  f = f*gain;
  % Apply analysis frame
  c = blockana(F, f); 
  % Apply thresholding
  c = thresh(c,thres,'soft');
  % Apply synthesis frame
  fhat = real(blocksyn(Fdual, c, bufLen));
  % Play the block
  blockplay(fhat);
end
% Close the control panel
p.close();
