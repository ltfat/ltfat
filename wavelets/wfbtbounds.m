function [AF,BF]=wfbtbounds(wt,L)
%WFBTBOUNDS Frame bounds of WFBT
%   Usage: fcond=wfbtbounds(wt,L);
%          [A,B]=wfbtbounds(wt,L);
%          [...]=wfbtbounds(wt);
%
%   `wfbtbounds(wt,L)` calculates the ratio $B/A$ of the frame bounds
%   of the filterbank tree specified by *wt* for a system of length
%   *L*. The ratio is a measure of the stability of the system.
%
%   `wfbtbounds({w,J,'dwt'},L)` calculates the ratio $B/A$ of the frame
%   bounds of the DWT (|FWT|) filterbank specified by *w* and *J* for a
%   system of length *L*.
%
%   `wfbtbounds(wt)` does the same thing, but *L* is assumed to be the
%   next compatible length bigger than the longest filter in the identical
%   filterbank.
%
%   `[A,B]=wfbtbounds(...)` returns the lower and upper frame bounds
%   explicitly.
%
%   See |wfbt| for explanation of parameter *wt* and |fwt| for explanation
%   of parameters *w* and *J*.
%
%   See also: wfbt, fwt, filterbankbounds

% AUTHOR: Zdenek Prusa

complainif_notenoughargs(nargin,1,'WFBTBOUNDS');

wt = wfbtinit({'strict',wt},'nat');

if nargin<2 
    L = [];
end;

if ~isempty(L)
   if L~=wfbtlength(L,wt)
       error(['%s: Specified length L is incompatible with the length of ' ...
              'the time shifts.'],upper(mfilename));
   end;
end

% Do the equivalent filterbank using multirate identity property
[gmultid,amultid] = wfbt2filterbank(wt);

if isempty(L)
   L = wfbtlength(max(cellfun(@(gEl) numel(gEl.h),gmultid)),wt);  
end

% Do the equivalent uniform filterbank
[gu,au] = nonu2ufilterbank(gmultid,amultid);

if nargout<2
   AF = filterbankbounds(gu,au,L);
elseif nargout == 2
   [AF, BF] = filterbankbounds(gu,au,L);
end
