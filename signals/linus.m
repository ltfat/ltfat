function [s,fs]=linus()
%LINUS  Load the 'linus' test signal
%   Usage:  s=linus;
%
%   `linus` loads the 'linus' signal. It is a recording of Linus Thorvalds
%   pronouncing the words "Hello. My name is Linus Thorvalds, and I
%   pronounce Linux as Linux".
%
%   The signal is 41461 samples long and is sampled at 8 kHz.
%
%   `[sig,fs]=linus` additionally returns the sampling frequency *fs*.
%
%   See `<http://www.paul.sladen.org/pronunciation/>`_.
  
%   AUTHOR : Peter L. Søndergaard
%   TESTING: TEST_SIGNALS
%   REFERENCE: OK

if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

s=wavload([f,'.wav']);
fs=8000;
