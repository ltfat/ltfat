function blockdone(varargin)
%BLOCKDONE  Destroy the current blockstream
%   Usage: blockdone();
%
%   `blockdone()` closes the current blockstream. The function resets
%   the playrec tool and clear all buffers in block_interface.
%
%   `blockdone(p1,p2,...)` in addition tries to call close methods on
%   all input arguments which are JAVA objects (which are passed by reference).
%   
%
%   See also: block

% TO DO: Process additional zeros to compensate for the delay 

block_interface('clearAll');
if playrec('isInitialised')
   playrec('reset');
end
clear playrec;

for ii=1:numel(varargin)
   p = varargin{ii};
   if isjava(p)
      try
         javaMethod('close',p);
      catch
         warning(sprintf('%s: Object %i does not have a close method.',...
         upper(mfilename),ii));
      end
   elseif isstruct(p) && isfield(p,'destructor') &&...
          isa(p.destructor,'function_handle')
      p.destructor();
   end
end
