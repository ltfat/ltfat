function [s,fs]=otoclick()
%OTOCLICK  Load the 'otoclick' test signal
%   Usage:  s=otoclick;
%
%   `otoclick` loads the 'otoclick' signal. The signal is a click-evoked
%   otoacoustic emission. It consists of two clear clicks followed by a
%   ringing. The ringing is the actual otoacoustic emission.
%
%   `[sig,fs]=otoclick` additionally returns the sampling frequency fs.
%
%   It was measured by Sarah Verhulst at CAHR (Centre of Applied Hearing
%   Research) at Department of Eletrical Engineering, Technical University
%   of Denmark
%
%   The signal is 2210 samples long and sampled at 44.1 kHz.

%   AUTHOR : Peter L. Søndergaard
%   TESTING: TEST_SIGNALS
%   REFERENCE: OK
  
if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

s=load('-ascii',[f,'.asc']);
fs=44100;

