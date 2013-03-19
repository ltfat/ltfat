function [y,n] = erbspacebw(flow,fhigh,varargin)
%ERBSPACEBW  Erbscale points specified by bandwidth
%   Usage: y=erbspacebw(flow,fhigh,bw,hitme);
%          y=erbspacebw(flow,fhigh,bw);
%          y=erbspacebw(flow,fhigh);
%
%   This is a wrapper around |audspacebw| that selects the erb-scale. Please
%   see the help on |audspacebw| for more information.
%
%   See also: audspacebw, freqtoaud
  
%   AUTHOR : Peter L. Søndergaard
  
[y,n] = audspacebw(flow,fhigh,varargin{:},'erb');
