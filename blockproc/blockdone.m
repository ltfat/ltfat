function blockdone(B)
%BLOCKDONE  Destroy the block object
%   Usage: blockdone(B);
%
%   `blockdone(B)` closes the block interface.
%
%   See also: block

% TO DO: Process additional zeros to compensate for the delay 

block_interface('free');