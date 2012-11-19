function freq = erbtofreq(erb);
%ERBTOFREQ  Converts erb units to frequency (Hz)
%   Usage: freq = erbtofreq(erb);
%  
%   This is a wrapper around |audtofreq|_ that selects the erb-scale. Please
%   see the help on |audtofreq|_ for more information.
%
%   See also: audtofreq, freqtoaud

%   AUTHOR: Peter L. Søndergaard
  
freq = audtofreq(erb,'erb');
