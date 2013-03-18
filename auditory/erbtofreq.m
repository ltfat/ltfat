function freq = erbtofreq(erb);
%ERBTOFREQ  Converts erb units to frequency (Hz)
%   Usage: freq = erbtofreq(erb);
%  
%   This is a wrapper around |audtofreq| that selects the erb-scale. Please
%   see the help on |audtofreq| for more information.
%
%   The following figure shows the corresponding frequencies for erb
%   values up to 31:::
%
%     erbs=0:31;
%     freqs=erbtofreq(erbs);
%     plot(erbs,freqs);
%     xlabel('Frequency / erb');
%     ylabel('Frequency / Hz');
%  
%   See also: audtofreq, freqtoaud

%   AUTHOR: Peter L. Søndergaard
  
freq = audtofreq(erb,'erb');
