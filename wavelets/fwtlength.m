function L=fwtlength(Ls,J);
%FWTLENGTH  FWT length from signal
%   Usage: L=fwtlength(Ls,J)
%
%   `L=fwtlength(Ls,J)` returns the length of the signal for non-expansive
%   discrete wavelet transform with periodic extensions
%
%   See also: dgt



if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;


L=2^J*ceil(Ls/2^J);
