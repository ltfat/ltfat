function bw = audfiltbw(fc)
%AUDFILTBW  Bandwidth of auditory filter
%   Usage: bw = audfiltbw(fc)
%
%   `audfiltbw(fc)` returns the equivalent rectangular bandwidth of the
%   auditory filter at center frequency *fc*. The function uses the
%   relation
%
%   ..  bw = 24.7 + fc/9.265
%
%   ..  math::  bw = 24.7 + \frac{fc}{9.265}
%     
%   as estimated in Glasberg and Moore (1990)
%
%   See also: freqtoerb, erbspace
%
%   References: glasberg1990daf
  
%   AUTHOR : Peter L. Søndergaard
  
% ------ Checking of input parameters ---------
  
error(nargchk(1,1,nargin));

if ~isnumeric(fc) || any(fc(:)<0)
  error('AUDFILTBW: fc must be non-negative.');
end;

% ------ Computation --------------------------

% FIXME: What is the upper frequency for which the estimation is valid?

bw = 24.7 + fc/9.265;

