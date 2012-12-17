function d=wilframediag(g,M,L,varargin)
%WILFRAMEDIAG  Diagonal of Wilson and WMDCT frame operator
%   Usage:  d=wilframediag(g,M,L);
%
%   Input parameters:
%         g     : Window function.
%         M     : Number of channels.
%         L     : Length of transform to do.
%   Output parameters:
%         d     : Diagonal stored as a column vector
%
%   `wilframediag(g,M,L)` computes the diagonal of the Wilson or WMDCT frame
%   operator with respect to the window *g* and number of channels *M*. The
%   diagonal is stored a as column vector of length *L*.
%
%   The diagonal of the frame operator can for instance be used as a
%   preconditioner.
%
%   See also: dwilt, wmdct, gabframediag

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

d=gabframediag(g,M,2*M,L)/2;


