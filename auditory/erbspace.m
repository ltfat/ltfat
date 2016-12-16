function [y,bw] = erbspace(fmin,fmax,n)
%ERBSPACE  Equidistantly spaced points on erbscale
%   Usage: y=erbspace(fmin,fmax,n);
%
%   This is a wrapper around |audspace| that selects the erb-scale. Please
%   see the help on |audspace| for more information.
%
%   See also: audspace, freqtoaud

%   AUTHOR : Peter L. Søndergaard
  
[y,bw] = audspace(fmin,fmax,n,'erb');

