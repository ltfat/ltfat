function freq = erbtofreq(erb);
%ERBTOFREQ  Converts erb units to frequency (Hz)
%   Usage: freq = erbtofreq(erb);
%  
%   This is a wrapper around audtofreq that selects the erb-scale. Please
%   see the help on AUDTOFREQ for more information.
%
%   See also: audtofreq, freqtoaud

%   AUTHOR: Peter L. Soendergaard
  
freq = audtofreq(erb,'erb');


