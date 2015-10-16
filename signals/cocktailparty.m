function [s,fs]=cocktailparty()
%COCKTAILPARTY  Load the 'cocktailparty' test signal
%   Usage:  s=cocktailparty;
%
%   `cocktailparty` loads the 'cocktailparty' signal. It is a recording of a
%   male native English speaker pronouncing the sentence "The cocktail party
%   effect refers to the ability to focus on a single talker among a mixture
%   of conversations in background noises".
%
%   `[sig,fs]=cocktailparty` additionally returns the sampling frequency *fs*.
%
%   The signal is 363200 samples long and recorded at 44.1 kHz in an
%   anechoic environment.

%   AUTHOR : James harte and Peter L. Søndergaard
%   TESTING: TEST_SIGNALS
%   REFERENCE: OK

if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

s=wavload([f,'.wav']);
fs=44100;
