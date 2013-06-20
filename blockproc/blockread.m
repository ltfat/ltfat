function [f,valid] = blockread(L)
%BLOCKREAD Read one block from input
%   Usage: blockread(L)
%       
%   Input parameters:
%      L    : Number of samples.
%   Output parameters:
%      f     : Samples.
%      valid : Input data valid flag.
%
%   Function also control the playback, so it does not have to rely on
%   whether user called |blockplay|.
% 
%   Block streaming uses several buffers to compensate for the processing
%   delay variation. 

persistent Lwav;
persistent clearStr;
persistent readTime;
%global delayLog;
%global delayLog2;

if nargin<1
   L = block_interface('getDefaultBufLen'); 
end


if block_interface('getDispLoad') 
   if block_interface('getPageNo')>0
      procTime = toc;
      res = 2;


   fs= playrec('getSampleRate');
   %delayLog = [delayLog, procTime];
   load = floor(100*(procTime+readTime)/(L/fs));
   msg = sprintf(['Load : |',repmat('*',1,ceil(min([load,100])/res)),repmat(' ',1,floor((100-min([load,100]))/res)),'| \n']);
   droppedStr = sprintf('Dropped samples: %i\n',playrec('getSkippedSampleCount'));
   fprintf([clearStr,msg,droppedStr]);
   clearStr = repmat(sprintf('\b'), 1, length(msg)+length(droppedStr));
   block_interface('setSkipped',playrec('getSkippedSampleCount'));
   if playrec('getSkippedSampleCount') > block_interface('getSkipped')
      block_interface('setSkipped',playrec('getSkippedSampleCount'));
   end

   else
      clearStr = '';
      %delayLog = [];
      %delayLog2 = [0];
      procTime = 0;
   end
   tic;
end

valid = 1;
source = block_interface('getSource');
pos = block_interface('getPos') +1; % convert to the Matlab indexing
block_interface('incPageNo');
pageNo = block_interface('getPageNo');
classid = block_interface('getClassId');

% Update sample counter
block_interface('setPos',pos+L-1); % convert back the Matlab indexing

%%%
%% REC, source is a mic/aux, no loopback
%

if strcmp(source,'rec')
   recChanList = block_interface('getRecChanList');
   
   readTime = toc;
   % Issue reading buffers up to max
   while block_interface('getEnqBufCount') <= block_interface('getBufCount')
      block_interface('pushPage', playrec('rec', L, recChanList));
   end
   pageList = block_interface('getPageList');
   % Block until the first page is loaded
   while(playrec('isFinished', pageList(1)) == 0)
   end
   % Read the data. Cast to the specified type
   f = cast(playrec('getRec',pageList(1)),classid);
   % Delete page
   playrec('delPage', pageList(1));
   % Throw away the page id
   block_interface('popPage');
   
%%%   
%% PLAYREC, source is a mic, loopback to an output
%
elseif strcmp(source,'playrec')
   recChanList = block_interface('getRecChanList');
   if pageNo<=1
      blockplay(zeros(L,numel(recChanList),classid));
   end
   
   % Enqueue already processed
   fhat = block_interface('getEnqueuedToPlay');
   if isempty(fhat)
      fhat = zeros(L,numel(recChanList),classid);
   end
   chanList = block_interface('getPlayChanList');

   % Copy input channel to all output chanels.
   fhat = repmat(fhat,1,numel(chanList));
   % Play and record
   block_interface('pushPage',playrec('playrec', fhat, chanList, -1, recChanList));

  readTime = toc;
   pageList = block_interface('getPageList');
   % Playback is block_interface('getBufCount') behind the input
   if block_interface('getPageNo') <= block_interface('getBufCount')
      f = zeros(L,numel(recChanList),classid);
   else
      % Block until the first page is loaded
      while(playrec('isFinished', pageList(1)) == 0)
      end
      % Read the data
      f = cast(playrec('getRec',pageList(1)),classid);
      playrec('delPage', pageList(1));
      % Throw away the page id
      block_interface('popPage');
   end

%%%   
%% PLAY: Source is a *.wav file
%
elseif strcmp(source(end-3:end),'.wav')
   % Get play channel list (could be chached) 
   chanList = block_interface('getPlayChanList');
   % Get already processed (from blockplay)
   fhat = block_interface('getEnqueuedToPlay');

   % Create something if blockplay was not called
   if isempty(fhat)
      fhat = zeros(L,numel(chanList),classid);
   end

   % Broadcast single input channel to all output chanels.
   if size(fhat,2)==1
      fhat = repmat(fhat,1,numel(chanList));
   end

   % Number of wav samples (is chached, since it is disk read operation)
   if pageNo<=1
      Lwav = wavread(source,'size'); 
   end

   if pos>Lwav(1)
      % Produce zeros when outside of the wav
      f = zeros(L,Lwav(2),classid);
      valid = 0;
   else
      % Determine valid samples
      endSample = min(pos + L - 1, Lwav(1));
      f = cast(wavread(source,[pos, endSample]),block_interface('getClassId')); 
      % Pad with zeros if some samples are missing
      if (pos + L - 1) > Lwav(1)
         ftmp = zeros(L,Lwav(2),classid);
         ftmp(1:size(f,1),:) = f;
         f = ftmp;
      end
   end

   % playrec('play',... - enques fhat to be played
   % block_interface('pushPage', - stores page number in an inner FIFO
   % queue
   block_interface('pushPage', playrec('play', fhat, chanList));

   readTime = toc;
   % If enough buffers are enqued, block the execution here until the 
   % first one is finished.
   if block_interface('getEnqBufCount') > block_interface('getBufCount')
      pageId = block_interface('popPage');
      % "Aggresive" chceking if page was played.
      % Another (supposedly slower) option is:
      % playrec('block',pageId);
      while(playrec('isFinished', pageId) == 0), end;
   end
end

if pageNo<=1
   playrec('resetSkippedSampleCount');
end

if block_interface('getDispLoad')
   tic;
end

