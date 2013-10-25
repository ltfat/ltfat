function blockdone()
%BLOCKDONE  Destroy the block object
%   Usage: blockdone(B);
%
%   `blockdone()` closes the current block interface.
%
%   See also: block

% TO DO: Process additional zeros to compensate for the delay 

block_interface('reset');
playrec('reset');