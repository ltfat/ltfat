function [y,bw] = erbspace(flow,fhigh,n)
%ERBSPACE  Equidistantly spaced points on erbscale
%   Usage: y=erbspace(flow,fhigh,n);
%
%   This is a wrapper around AUDSPACE that selects the erb-scale. Please
%   see the help on AUDSPACE for more information.
%
%   See also: audspace, freqtoaud

%   AUTHOR : Peter L. Soendergaard
  
[y,bw] = audspace(flow,fhigh,n,'erb');
  
