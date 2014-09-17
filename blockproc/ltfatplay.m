function ltfatplay(source,varargin)
%LTFATPLAY Play data samples or a wav file
%   Usage: ltfatplay('file.wav')
%          ltfatplay(data,'fs',fs)
%          ltfatplay(...,'devid',devid)
%
%
%   `ltfatplay('file.wav')` plays a wav file using the default sound device.
%
%   `ltfatplay('file.wav','devid',devid)` plays a wav file using the sound
%   device with id `devid`. A list of available devices can be obtained by 
%   |blockdevices|.
%
%   `ltfatplay(data,'fs',fs,...)` works entirely similar, but `data` is
%   expected to be a vector of length $L$ or a $L\times W$ matrix with
%   columns as individual channels and *fs* to be a sampling rate to be used.
%   When no sampling rate is specified, 44.1 kHz is used.
%
%   In addition, individual channels of the output sound device can be
%   selected by using an additional key-value pair
%
%   `'playch',playch`
%      A vector of channel indexes starting at 1.
%
%   This function has the advantage over `sound` and `soundsc` that one can 
%   directly specify output device chosen from |blockdevices|. Similar
%   behavior can be achieved using `audioplayer` and `audiodevinfo` but
%   only in Matlab. Audioplayer is not yet supported in Octave.
%   

%   Author: Zdenek Prusa

% Initialize block stream
block(source,varargin{:});

sourceHandle = block_interface('getSource');

% Allow playing data and wav files only
if ~isa(sourceHandle,'function_handle')
    error('%s: Specified source cannot be played.',upper(mfilename));
end

% Make the vector safe
Ls = block_interface('getLs');
Lssafe = max([256,Ls(1)]);
f = postpad(sourceHandle(1,Ls(1)),Lssafe);

% If one channel is used, broadcast it to all output channels
chanList = block_interface('getPlayChanList');
if size(f,2)==1
    f = repmat(f,1,numel(chanList));
end

% Finally play it at once
playrec('play',f,chanList);



