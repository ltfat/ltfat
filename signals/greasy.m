function [s,fs]=greasy()
%GREASY  Load the 'greasy' test signal.
%   Usage:  s=greasy;
%
%   GREASY loads the 'greasy' signal. It is a recording of a woman
%   pronouncing the word "greasy".
%
%   The signal is 5880 samples long and recorded at 16 kHz with around 11
%   bits of effective quantization.
%
%   [sig,fs]=GREASY additionally returns the sampling frequency fs.
%
%   The signal has been scaled to not produce any clipping when
%   played. To get integer values use round(greasy*2048).
%
%   The signal was obtained from Wavelab:
%     http://www-stat.stanford.edu/~wavelab/
%
%R  mazh93

%   AUTHOR : Peter Soendergaard
%   TESTING: TEST_SIGNALS
%   REFERENCE: OK
  
if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

s = wavread([f,'.wav']);
fs = 16000;
