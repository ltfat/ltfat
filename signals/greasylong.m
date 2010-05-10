function s=greasylong()
%GREASYLONG  Load the 'greasylong' test signal.
%   Usage:  s=greasylong;
%
%   GREASYLONG loads the 'greasylong' signal. It is a recording of a woman
%   pronouncing the sentence 'She had your dark suit in greasy wash-water
%   all year'. The sentence is the very first one in the TIMIT database.
%
%   The signal is 20768 samples long and recorded at 8 khz.
%
%   The signal was obtained from the CAAR publications page:
%     http://www.isr.umd.edu/CAAR/pubs.html

%   AUTHOR : Peter Soendergaard
%   TESTING: TEST_SIGNALS
%   REFERENCE: OK

if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

s=wavread([f,'.wav']);


