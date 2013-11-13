function [AF,BF]=fwtbounds(w,J,L)
%FWTBOUNDS Frame bounds of DWT
%   Usage: fcond=fwtbounds(w,J,L);
%          [A,B]=fwtbounds(w,J,L);
%
%   `fwtbounds(w,a,L)` calculates the ratio $B/A$ of the frame bounds
%   of the filterbank specified by *w* and *J* for a system of length
%   *L*. The ratio is a measure of the stability of the system. 
%
%   `[A,B]=fwtbounds(w,J,L)` returns the lower and upper frame bounds
%   explicitly. 
%
%   See |fwt| for explanation of parameters *w* and *J*.
%
%   See also: fwt, filterbankbounds


if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

% There seem to be two possibilites:
% 1) Find the frame bounds of the equaivalent uniform filterbank. The
%    equivalent non-uniform filterbank can be created using fwt2filterbank,
%    non-uniform to uniform filterbank transform is described in the book:
%    Beyond Wavelets, chapter 10, Nonuniform filter banks: New results and
%    open problems by Akkarakaran and Vaidyanathan.
% 2) Use the recursive polyphase matrix multiplication algorithm from 
%    Stanhill, D.; Zeevi, Y. Y.: Frame analysis of wavelet-type filter banks 

%ad 1) (shotcut trough wfbtbounds)

if nargout<2
   AF = wfbtbounds({w,J,'dwt'},L);
elseif nargout == 2
   [AF, BF] = wfbtbounds({w,J,'dwt'},L);
end
