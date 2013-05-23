function blockplay(f)
%BLOCKPLAY Plays block
%
%
%

source = block_interface('getSource');

if strcmp(source,'rec')
   error('%s: Blocks cannot be played in the rec only mode.',upper(mfilename));
end

block_interface('enqueueToPlay',f);

% source = block_interface('getSource');
% chanList = block_interface('getPlayChanList');
% 
% if strcmp(source,'playrec')
%    recChanList = block_interface('getRecChanList');
%    block_interface('pushPage',playrec('playrec', f, chanList, -1, recChanList));
%    % Blocking is handled by the input
% else
%    block_interface('pushPage', playrec('play', f, chanList));
% 
%    % If enough buffers are enqued, block the execution here until the first one is finished
%      if block_interface('getEnqBufCount') > block_interface('getBufCount')
%         pageId = block_interface('popPage');
%         while(playrec('isFinished', pageId) == 0)
%         end
%      end
% end