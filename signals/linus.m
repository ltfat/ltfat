function s=linus()
%LINUS  Load the 'linus' test signal.
%   Usage:  s=linus;
%
%   LINUS loads the 'linus' signal. It is a recording
%   of Linus Thorvalds pronouncing the words "Hello. My name is Linus
%   Thorvalds, and I pronounce Linux as Linux".
%
%   The signal is 41461 samples long and is sampled at 8kHz.
%
%   See http://www.paul.sladen.org/pronunciation/
  
%   AUTHOR : Peter Soendergaard
%   TESTING: TEST_SIGNALS
%   REFERENCE: OK

if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

s=wavread([f,'.wav']);

