function [AF,BF]=ufwtbounds(w,J,L,varargin)
%UFWTBOUNDS Frame bounds of Undecimated DWT
%   Usage: fcond=ufwtbounds(w,J,L);
%          [A,B]=ufwtbounds(w,J,L);
%
%   `ufwtbounds(w,a,L)` calculates the ratio $B/A$ of the frame bounds
%   of the filterbank specified by *w* and *J* for a system of length
%   *L*. The ratio is a measure of the stability of the system.
%
%   `[A,B]=ufwtbounds(w,J,L)` returns the lower and upper frame bounds
%   explicitly.
%
%   See |ufwt| for explanation of parameters *w* and *J*.
%
%   See also: ufwt, filterbankbounds


complainif_notenoughargs(nargin,3,'UFWTBOUNDS');

if nargout<2
   AF = uwfbtbounds({w,J,'dwt'},L,varargin{:});
elseif nargout == 2
   [AF, BF] = uwfbtbounds({w,J,'dwt'},L,varargin{:});
end
