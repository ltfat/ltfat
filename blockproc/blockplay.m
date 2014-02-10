function blockplay(f)
%BLOCKPLAY Schedules block to be played
%   Usage: blockplay(L)
%       
%   Input parameters:
%      f    : Samples.
%
%   Function schedules samples in *f* to be played. Since playrec handles
%   playing and recording in a single command, the actual relay of samples
%   to playrec is done in the next call of |blockread|.


source = block_interface('getSource');

if strcmp(source,'rec')
   error('%s: Blocks cannot be played in the rec only mode.',upper(mfilename));
end

block_interface('setToPlay',f);
