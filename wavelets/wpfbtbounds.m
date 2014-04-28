function [AF,BF]=wpfbtbounds(wt,L,varargin)
%WPFBTBOUNDS Frame bounds of WPFBT
%   Usage: fcond=wpfbtbounds(wt,L);
%          [A,B]=wpfbtbounds(wt,L);
%
%   `wpfbtbounds(wt,L)` calculates the ratio $B/A$ of the frame bounds
%   of the wavelet packet filterbank specified by *wt* for a system of length
%   *L*. The ratio is a measure of the stability of the system.
%
%   `[A,B]=wpfbtbounds(wt,L)` returns the lower and upper frame bounds
%   explicitly.
%
%   See |wfbt| for explanation of parameter *wt*.
%
%   See also: wpfbt, filterbankbounds


complainif_notenoughargs(nargin,2,'WPFBTBOUNDS');

definput.flags.interscaling = {'intsqrt', 'intscale', 'intnoscale'};
[flags]=ltfatarghelper({},definput,varargin);

wt = wfbtinit({'strict',wt},'nat');

if L~=wfbtlength(L,wt)
    error(['%s: Specified length L is incompatible with the length of ' ...
           'the time shifts.'],upper(mfilename));
end;


for ii=1:numel(wt.nodes)
   a = wt.nodes{ii}.a;
   assert(all(a==a(1)),sprintf(['%s: One of the basic wavelet ',...
                                'filterbanks is not uniform.'],...
                                upper(mfilename)));
end

% Scale the intermediate filters
wt = comp_wpfbtscale(wt,flags.interscaling);

% Do the equivalent filterbank using multirate identity property
[gmultid,amultid] = wpfbt2filterbank(wt);

% Do the equivalent uniform filterbank
[gu,au] = nonu2ufilterbank(gmultid,amultid);

if nargout<2
   AF = filterbankbounds(gu,au,L);
elseif nargout == 2
   [AF, BF] = filterbankbounds(gu,au,L);
end
