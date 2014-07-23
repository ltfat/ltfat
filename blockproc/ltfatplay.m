function ltfatplay(source,varargin)
%LTFATPLAY Play data samples or a wav file
%
%
%
%   `ltfatplay(source)` plays data or a wav file specified as `source`.
%   This function has advantage over sound and soundsc that one can
%   directly specify output device chosen from |blockdevices|. Similar
%   behavior can be achieved using audioplayer and audiodevinfo in Matlab
%   only. Audioplayer is not yet supported in Octave.
%   



block(source,varargin{:});

sourceHandle = block_interface('getSource');

if ~isa(sourceHandle,'function_handle')
    error('%s: Specified source cannot be played.',upper(mfilename));
end

Ls = block_interface('getLs');
Lssafe = max([256,Ls(1)]);
f = postpad(sourceHandle(1,Ls(1)),Lssafe);

chanList = block_interface('getPlayChanList');
if size(f,2)==1
    f = repmat(f,1,numel(chanList));
end

playrec('play',f,chanList);



