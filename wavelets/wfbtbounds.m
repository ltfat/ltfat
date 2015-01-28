function [AF,BF]=wfbtbounds(wt,varargin)
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
%   The function supports the following flag groups:
%
%   `'scaling_notset'`(default),`'noscale'`,`'scale'`,`'sqrt'`
%     Support for scaling flags as described in |uwfbt|. By default,
%     the bounds are computed for |wfbt|, passing any of the non-default
%     flags results in framebounds for |uwfbt|.
%
%   See also: wfbt, fwt, filterbankbounds

% AUTHOR: Zdenek Prusa

complainif_notenoughargs(nargin,1,'WFBTBOUNDS');

% Frequency or natural ordering should not make a difference
wt = wfbtinit({'strict',wt},'nat');

definput.keyvals.L = [];
definput.import = {'uwfbtcommon'};
definput.importdefaults = {'scaling_notset'};
[flags,~,L]=ltfatarghelper({'L'},definput,varargin);

if ~isempty(L) && flags.do_scaling_notset
   if L~=wfbtlength(L,wt)
       error(['%s: Specified length L is incompatible with the length of ' ...
              'the time shifts.'],upper(mfilename));
   end;
end

% Do the equivalent filterbank using multirate identity property
[gmultid,amultid] = wfbt2filterbank(wt,flags.scaling);

if isempty(L)
   L = wfbtlength(max(cellfun(@(gEl) numel(gEl.h),gmultid)),wt);  
end

% Do the equivalent uniform filterbank
if any(amultid~=amultid(1))
   [gu,au] = nonu2ufilterbank(gmultid,amultid);
else
   [gu,au] = deal(gmultid,amultid);
end

if nargout<2
   AF = filterbankbounds(gu,au,L);
elseif nargout == 2
   [AF, BF] = filterbankbounds(gu,au,L);
end
