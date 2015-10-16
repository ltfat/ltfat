function [s,fs]=greasy()
%GREASY  Load the 'greasy' test signal
%   Usage:  s=greasy;
%
%   `greasy` loads the 'greasy' signal. It is a recording of a woman
%   pronouncing the word "greasy".
%
%   The signal is 5880 samples long and recorded at 16 kHz with around 11
%   bits of effective quantization.
%
%   `[sig,fs]=greasy` additionally returns the sampling frequency *fs*.
%
%   The signal has been scaled to not produce any clipping when
%   played. To get integer values use `round(greasy*2048)`.
%
%   The signal was obtained from Wavelab:
%   `<http://www-stat.stanford.edu/~wavelab/>`_, it is a part of the first
%   sentence of the TIMIT speech corpus "She had your dark suit in greasy
%   wash water all year":
%   `<http://www.ldc.upenn.edu/Catalog/CatalogEntry.jsp?catalogId=LDC93S1>`_.
%
%   Examples:
%   ---------
%
%   Plot of 'greasy' in the time-domain:::
%
%     plot((1:5880)/16000,greasy);
%     xlabel('Time (seconds)');
%     ylabel('Amplitude');
%
%   Plot of 'greasy' in the frequency-domain:::
%
%     plotfftreal(fftreal(greasy),16000,90);
%
%   Plot of 'greasy' in the time-frequency-domain:::
%
%     sgram(greasy,16000,90);
%
%   References: mazh93

%   AUTHOR : Peter L. Søndergaard
%   TESTING: TEST_SIGNALS
  
if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

s = wavload([f,'.wav']);
fs = 16000;

