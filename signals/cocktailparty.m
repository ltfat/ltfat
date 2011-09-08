function s=cocktailparty()
%COCKTAILPARTY  Load the 'cocktailparty' test signal.
%   Usage:  s=cocktailparty;
%
%   COCKTAILPARTY loads the 'cocktailparty' signal. It is a recording of a
%   male native English speaker pronouncing the sentence 'The cocktail party
%   effect refers to the ability to focus on a single talker among a mixture
%   of conversations in background noises'.
%
%   The signal is 363200 samples long and recorded at 44100 khz in an
%   anechoic environment.

%   AUTHOR : James harte and Peter Soendergaard
%   TESTING: TEST_SIGNALS
%   REFERENCE: OK

if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

s=wavread([f,'.wav']);


