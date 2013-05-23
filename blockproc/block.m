function [fs,classid] = block(source,varargin)
%BLOCKINIT Initialize block stream
%   Usage: block(source)
%          block(source,'fs',fs,'nbuf',nbuf,'fmt',fmt)  
%
%   Input parameters:
%      source    : Block stream input.
%      fs        : Sampling rate.
%      nbuf      : Number of buffers.
%      fmt       : Data type
%   Output parameters:
%      fs        : Sampling rate.
%      classid   : Data type.
%
%   `block(source)` initializes block data stream from `source` which
%   can be (the letter-case is ignored for strings):
%
%      `'file.wav'`    : name of a wav file
%
%      `'dialog'`      : shows the file dialog to choose a wav file.
%
%      `'rec'`         : input is taken from a microfone/auxilary input
%
%      `'playrec'`     : loopbacks the input to the output
%
%      `data`          : input data as columns of a matrix for each input 
%                        channel
%
%   `block` accepts the following optional flags and key-value pairs:
%
%   Optional key-value pairs:
%
%      `'devid',dev`   : Whenever more input/output devices are present in
%                        your system, `'devid'` can be used to specify one.
%                        For the `'playrec'` option the `devId` should be a
%                        two element vector [playDevid, recDevid]. List
%                        of the installed devices ant their ids can be 
%                        obtained by |block_devices|.
%      
%      `'playch',playch`: If device supports more output channels, `'playch'`
%                         can be used to specify which should be used. E.g.
%                         for two channel device, [1,2] is used to specify 
%                         both channels. 
%
%      `'recch',recch`  : If device supports more input channels, `'recch'`
%                         can be used to specify which should be used.
%
%      `'outfile','file.wav'`: Processed sound data is stored in a new wav
%                              file.


definput.keyvals.devid=[];
definput.keyvals.nbuf=[3];
definput.keyvals.fs=44100;
definput.keyvals.playch=[];
definput.keyvals.recch=[];
definput.keyvals.outfile=[];
definput.flags.fmt={'double','single'};
[flags,kv]=ltfatarghelper({},definput,varargin);

playChannels = 0;
recChannels = 0;
play = 0;
record = 0;
paBufsize = 256;

if ~isempty(kv.outfile)
   error('%s: TO DO: Writing to the output wav file is not supported yet.',upper(mfielname));
end

if ischar(source)
   if(strcmpi(source,'rec'))
      recChannels = 1;
      record = 1;
   elseif strcmpi(source,'playrec')
      playChannels = 2;
      recChannels = 1;
      record = 1;   
      play = 1;
   elseif strcmpi(source,'dialog')
      error('%s: TO DO: Open dialog for a sound file.',upper(mfilename));
   elseif(numel(source)>4)
      if(strcmp(source(end-3:end),'.wav'))
         if exist(source)~=2
            error('%s: File "%s" does not exist.',upper(mfilename),source);
         end
         [kv.L, kv.fs] = wavread(source, 'size');
         playChannels = 2;
         kv.L = kv.L(1);
         play = 1;
      else
         error('%s: "%s" is not valid wav filename.',upper(mfilename),source);
      end

   else
      error('%s: Unrecognized command "%s".',upper(mfilename),source);
   end
elseif(isnumeric(source))
   error('%s: TO DO: Accept samples as a vector or matrix.',upper(mfilename));
   kv.L = size(source,1);
   playChannels = 2;
   play = 1;
else
   error('%s: Unrecognized input.',upper(mfilename));
end

isPlayrecInit = 0;
try 
   isPlayrecInit = playrec('isInitialised');
catch
   err = lasterror;
   if ~isempty(strfind(err.message,'The specified module could not be found'))
      error('%s: playrec found but portaudio cannot be found.', upper(mfilename));
   end
    if ~isempty(strfind(err.message,'Undefined function'))
      error('%s: playrec could not be found.', upper(mfilename));
    end
   error('%s: Error loading playrec.',upper(mfilename));
end

if isPlayrecInit
   playrec('reset');
end

if isempty(kv.playch)
  kv.playch = 1:playChannels; 
end

if isempty(kv.recch)
  kv.recch = 1:recChannels; 
end

devs = playrec('getDevices');
if isempty(devs)
   error('%s: No sound devices available. portaudio lib is probably incorrectly built.',upper(mfilename));
end

block_interface('reset');
% Get all installed play devices
playDevIds = arrayfun(@(dEl) dEl.deviceID,devs(arrayfun(@(dEl) dEl.outputChans,devs)>0));
% Get all installed recording devices
recDevIds = arrayfun(@(dEl) dEl.deviceID,devs(arrayfun(@(dEl) dEl.inputChans,devs)>0));

if play && record
   if ~isempty(kv.devid)
      if(numel(kv.devid)~=2)
         error('%s: devid should be 2 element vector.',upper(mfilename));
      end
      if ~any(playDevIds==kv.devid(1))
         error('%s: There is no play device with id = %i.',upper(mfilename),kv.devid(1));
      end
      if ~any(recDevIds==kv.devid(2))
         error('%s: There is no rec device with id = %i.',upper(mfilename),kv.devid(2));
      end
   else
      kv.devid = [playDevIds(1), recDevIds(1)];
   end
   playrec('init', kv.fs, kv.devid(1), kv.devid(2), numel(kv.playch),numel(kv.recch),paBufsize);
   if numel(kv.recch) >1
      error('%s: Using more than one input channel.',upper(mfilename));
   end
   block_interface('setPlayChanList',kv.playch);
   block_interface('setRecChanList',kv.recch);
elseif play && ~record
   if ~isempty(kv.devid)
      if numel(kv.devid) >1
         error('%s: devid should be scalar.',upper(mfilename));
      end
      if ~any(playDevIds==kv.devid)
         error('%s: There is no play device with id = %i.',upper(mfilename),kv.devid);
      end
   else
      % Use the first (hopefully default) device
      kv.devid = playDevIds(1);
   end
   playrec('init', kv.fs, kv.devid, -1,numel(kv.playch),-1,paBufsize);
   block_interface('setPlayChanList',kv.playch);
   if(playrec('getPlayMaxChannel')<numel(kv.playch))
       error ('%s: Selected device does not support %d output channels.\n',upper(mfilename), max(chanList));
   end
elseif ~play && record
   if ~isempty(kv.devid)
      if ~any(recDevIds==kv.devid)
         error('%s: There is no rec device with id = %i.',upper(mfilename),kv.devid);
      end
   else
      % Use the first (hopefully default) device
      kv.devid = recDevIds(1);
   end
   playrec('init', kv.fs, -1, kv.devid,-1,numel(kv.recch),paBufsize);
   block_interface('setRecChanList',kv.recch);
else
   error('%s: Play or record should have been set.',upper(mfilename));
end

% From the playrec author:
% This slight delay is included because if a dialog box pops up during
% initialisation (eg MOTU telling you there are no MOTU devices
% attached) then without the delay Ctrl+C to stop playback sometimes
% doesn't work.
pause(0.1);   

if(~playrec('isInitialised'))
    error ('%s: Unable to initialise playrec correctly.',upper(mfilename));
end


if(playrec('pause'))
    %fprintf('Playrec was paused - clearing all previous pages and unpausing.\n');
    playrec('delPage');
    playrec('pause', 0);
end

% Reset skipped samples
playrec('resetSkippedSampleCount');

% Store data
block_interface('setSource',source);
block_interface('setBufCount',kv.nbuf);
block_interface('setClassId',flags.fmt);


% Return parameters
classid = flags.fmt;
fs = kv.fs;

   

      
   
   
   
   
   
   
   
   