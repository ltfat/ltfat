function erb = freqtoerb(freq);
%FREQTOERB  Converts frequencies (Hz) to erbs.
%   Usage: erb = freqtoerb(freq);
%
%   This is a wrapper around freqtoaud that selects the erb-scale. Please
%   see the help on FREQTOAUD for more information.
%
%   See also: freqtoaud
%
%   Demos: demo_audscales
  
%   AUTHOR: Peter L. Soendergaard

erb = freqtoaud(freq,'erb');