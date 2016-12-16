function [y,n] = erbspacebw(fmin,fmax,varargin)
%ERBSPACEBW  Erbscale points specified by bandwidth
%   Usage: y=erbspacebw(fmin,fmax,bw,hitme);
%          y=erbspacebw(fmin,fmax,bw);
%          y=erbspacebw(fmin,fmax);
%
%   This is a wrapper around |audspacebw| that selects the erb-scale. Please
%   see the help on |audspacebw| for more information.
%
%   See also: audspacebw, freqtoaud
  
%   AUTHOR : Peter L. Søndergaard
  
[y,n] = audspacebw(fmin,fmax,varargin{:},'erb');
