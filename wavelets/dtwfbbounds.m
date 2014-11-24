function [AF,BF]=dtwfbbounds(dualwt,L)
%DTWFBBOUNDS Frame bounds of DTWFB
%   Usage: fcond=dtwfbbounds(dualwt,L);
%          [A,B]=dtwfbbounds(dualwt,L);
%          [...]=dtwfbbounds(dualwt);
%
%   `dtwfbbounds(dualwt,L)` calculates the ratio $B/A$ of the frame bounds
%   of the dual-tree filterbank specified by *dualwt* for a system of 
%   length *L*. The ratio is a measure of the stability of the system.
%
%   `dtwfbbounds(dualwt)` does the same thing, but *L* is the next compatible 
%   length bigger than the longest filter in the identical filterbank.
%
%   `[A,B]=dtwfbbounds(...)` returns the lower and upper frame bounds
%   explicitly.
%
%   See |dtwfb| for explanation of parameter *dualwt*.
%
%   See also: dtwfb, filterbankbounds
%

% AUTHOR: Zdenek Prusa

complainif_notenoughargs(nargin,1,'DTWFBBOUNDS');

dualwt = dtwfbinit({'strict',dualwt},'nat');

if nargin<2 
    L = [];
end;

if ~isempty(L)
   if L~=wfbtlength(L,dualwt)
          error(['%s: Specified length L is incompatible with the length of ' ...
                 'the time shifts.'],upper(mfilename));
   end
end

% Do the equivalent filterbank using multirate identity property
[gmultid,amultid] = dtwfb2filterbank(dualwt,'complex');

if isempty(L)
   L = wfbtlength(max(cellfun(@(gEl) numel(gEl.h),gmultid)),dualwt);  
end


% Do the equivalent uniform filterbank
[gu,au] = nonu2ufilterbank(gmultid,amultid);

if nargout<2
   AF = filterbankbounds(gu,au,L);
elseif nargout == 2
   [AF, BF] = filterbankbounds(gu,au,L);
end

