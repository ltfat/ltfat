function c=batmask()
%BATMASK  Load a Gabor multiplier symbol for the 'bat' test signal
%   Usage:  c=batmask;
%
%   `batmask` loads a Gabor multiplier with a 0/1 symbol that masks out
%   the main contents of the 'bat' signal. The symbol fits a Gabor
%   multiplier with lattice given by $a=10$ and $M=40$.
%
%   The mask was created manually using a image processing program. The
%   mask is symmetric, such that the result will be real valued if the
%   multiplier is applied to a real valued signal using a real valued
%   window.
%
%   See also:  bat

%   AUTHOR : Peter L. Søndergaard
%   TESTING: TEST_BATMASK
%   REFERENCE: OK
  
if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

c=load('-ascii',[f,'.asc']);

