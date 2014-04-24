function [AF,BF]=uwfbtbounds(wt,L,varargin)
%WFBTBOUNDS Frame bounds of Undecimated WFBT
%   Usage: fcond=uwfbtbounds(wt,L);
%          [A,B]=uwfbtbounds(wt,L);
%
%   `uwfbtbounds(wt,L)` calculates the ratio $B/A$ of the frame bounds
%   of the undecimated filterbank specified by *wt* for a system of length
%   *L*. The ratio is a measure of the stability of the system.
%
%   `wfbtbounds({w,J,'dwt'},L)` calculates the ratio $B/A$ of the frame
%   bounds of the DWT (|FWT|) filterbank specified by *w* and *J* for a
%   system of length *L*.
%
%   `[A,B]=uwfbtbounds(...,L)` returns the lower and upper frame bounds
%   explicitly.
%
%   See |wfbt| for explanation of parameter *wt* and |fwt| for explanation
%   of parameters *w* and *J*.
%
%   See also: uwfbt, filterbankbounds


complainif_notenoughargs(nargin,2,'UWFBTBOUNDS');

definput.flags.scaling={'sqrt','scale','noscale'};
[flags]=ltfatarghelper({},definput,varargin);

wt = wfbtinit({'strict',wt},'nat');

for ii=1:numel(wt.nodes)
   a = wt.nodes{ii}.a;
   assert(all(a==a(1)),sprintf(['%s: One of the basic wavelet ',...
                                'filterbanks is not uniform.'],...
                                upper(mfilename)));
end

% Do the equivalent filterbank using multirate identity property
[gmultid, amultid] = wfbt2filterbank(wt);

% Scale filters
gmultid = comp_filterbankscale(gmultid, amultid, flags.scaling);

if nargout<2
   AF = filterbankbounds(gmultid,1,L);
elseif nargout == 2
   [AF, BF] = filterbankbounds(gmultid,1,L);
end
