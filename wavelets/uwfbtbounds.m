function [AF,BF]=uwfbtbounds(wt,varargin)
%UWFBTBOUNDS Frame bounds of Undecimated WFBT
%   Usage: fcond=uwfbtbounds(wt,L);
%          [A,B]=uwfbtbounds(wt,L);
%          [...]=uwfbtbounds(wt);
%
%   `uwfbtbounds(wt,L)` calculates the ratio $B/A$ of the frame bounds
%   of the undecimated filterbank specified by *wt* for a system of length
%   *L*. The ratio is a measure of the stability of the system.
%
%   `uwfbtbounds({w,J,'dwt'},L)` calculates the ratio $B/A$ of the frame
%   bounds of the undecimated DWT (|UFWT|) filterbank specified by *w* and
%   *J* for a system of length *L*.
%
%   `uwfbtbounds(wt)` does the same thing, but *L* is the length of the 
%   longest filter in the identical filterbank.
%
%   `[A,B]=uwfbtbounds(...)` returns the lower and upper frame bounds
%   explicitly.
%
%   See |wfbt| for explanation of parameter *wt* and |fwt| for explanation
%   of parameters *w* and *J*.
%
%   The function supports the following flags:
%
%   `'sqrt'`(default),`'noscale'`,`'scale'`
%       The filters in the filterbank tree are scaled to reflect the
%       behavior of |uwfbt| and |iuwfbt| with the same flags.  
%
%   See also: uwfbt, filterbankbounds


complainif_notenoughargs(nargin,1,'UWFBTBOUNDS');

definput.keyvals.L = [];
definput.flags.scaling={'sqrt','scale','noscale'};
[flags,~,L]=ltfatarghelper({'L'},definput,varargin);

wt = wfbtinit({'strict',wt},'nat');

for ii=1:numel(wt.nodes)
   a = wt.nodes{ii}.a;
   assert(all(a==a(1)),sprintf(['%s: One of the basic wavelet ',...
                                'filterbanks is not uniform.'],...
                                upper(mfilename)));
end

% Do the equivalent filterbank using multirate identity property
[gmultid, amultid] = wfbt2filterbank(wt);

if isempty(L)
   L = wfbtlength(max(cellfun(@(gEl) numel(gEl.h),gmultid)),wt);  
end

% Scale filters
gmultid = comp_filterbankscale(gmultid, amultid, flags.scaling);

if nargout<2
   AF = filterbankbounds(gmultid,1,L);
elseif nargout == 2
   [AF, BF] = filterbankbounds(gmultid,1,L);
end
