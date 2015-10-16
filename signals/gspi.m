function [s,fs]=gspi()
%GSPI  Load the 'glockenspiel' test signal
%
%   `gspi` loads the 'glockenspiel' signal. This is a recording of a simple
%   tune played on a glockenspiel. It is 262144 samples long, mono, recorded
%   at 44100 Hz using 16 bit quantization.
%   
%   `[sig,fs]=gspi` additionally returns the sampling frequency *fs*.
%
%   This signal, and other similar audio tests signals, can be found on
%   the EBU SQAM test signal CD `<http://tech.ebu.ch/publications/sqamcd>`_.
%
  
%   AUTHOR : Peter L. Søndergaard
%   TESTING: TEST_SIGNALS
%   REFERENCE: OK
  
if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

% Load audio signal
s = wavload([f,'.wav']);
fs = 44100;

