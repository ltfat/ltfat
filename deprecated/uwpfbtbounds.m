function [AF,BF]=uwpfbtbounds(wt,varargin)
%UWPFBTBOUNDS Frame bounds of Undecimated WPFBT
%   Usage: fcond=uwpfbtbounds(wt,L);
%          [A,B]=uwpfbtbounds(wt,L);
%
%   `uwpfbtbounds(wt,L)` calculates the ratio $B/A$ of the frame bounds
%   of the undecimated wavelet packet filterbank specified by *wt* for a 
%   system of length *L*. The ratio is a measure of the stability of the 
%   system. 
%
%   `[A,B]=uwfbtbounds(wt,L)` returns the lower and upper frame bounds
%   explicitly. 
%
%   See |wfbt| for explanation of parameter *wt*.
%
%   Additionally, the function accepts the following flags:
%
%   `'intsqrt'`(default),`'intnoscale'`, `'intscale'`
%       The filters in the filterbank tree are scaled to reflect the
%       behavior of |uwpfbt| and |iuwpfbt| with the same flags.
%
%   `'sqrt'`(default),`'noscale'`,`'scale'`
%       The filters in the filterbank tree are scaled to reflect the
%       behavior of |uwpfbt| and |iuwpfbt| with the same flags.  
%              
%
%   See also: uwpfbt, filterbankbounds

% AUTHOR: Zdenek Prusa

warning('UWPFBTBOUNDS is deprecated. Instead, use WPFBTBOUNDS with appropriate flags');

complainif_notenoughargs(nargin,1,'UWPFBTBOUNDS');

definput.keyvals.L = [];
definput.flags.scaling={'sqrt','scale','noscale'};
definput.flags.interscaling = {'intsqrt', 'intscale', 'intnoscale'};
[flags,~,L]=ltfatarghelper({'L'},definput,varargin);

if nargout<2
   AF = wpfbtbounds(wt,L,flags.scaling,flags.interscaling);
elseif nargout == 2
   [AF, BF] = wpfbtbounds(wt,L,flags.scaling,flags.interscaling);
end
