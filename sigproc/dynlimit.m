function xo=dynlimit(xi,dynrange,varargin);
%DYNLIMIT  Limit the dynamical range of the input
%   Usage: xo=dynlimit(xi,dynrange);
%
%   `dynlimit(xi,dynrange)` will threshold the input such that the
%   difference between the maximum and minumum value of *xi* is exactly
%   *dynrange*.
%
%   See also: thresh, largestr, largestn
  
xmax=max(xi(:));
xo=xi;
xo(xo<xmax-dynrange)=xmax-dynrange;
