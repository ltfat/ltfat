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
%   In case no audio output is expected (in the rec only mode), 
%   the function does nothing.


complainif_notenoughargs(nargin,1,'BLOCKPLAY');

source = block_interface('getSource');

if ( iscell(source) && strcmp(source{1},'rec')) || ...
   strcmp(source,'rec')
   % Do nothing in rec only mode.
   return; 
   % error('%s: Blocks cannot be played in the rec only mode.',upper(mfilename));
end

% Reformat f if necessary
f = comp_sigreshape_pre(f,'BLOCKPLAY',0);

block_interface('setToPlay',f);
