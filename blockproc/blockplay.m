function blockplay(f)
%BLOCKPLAY Plays block
%
%
%

source = block_interface('getSource');

if strcmp(source,'rec')
   error('%s: Blocks cannot be played in the rec only mode.',upper(mfilename));
end

block_interface('setToPlay',f);
