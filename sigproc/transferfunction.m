function H=transferfunction(g,L)
%TRANSFERFUNCTION  The transferfunction of a filter
%   Usage:  H=transferfunction(g,L);
%
%   `transferfunction(g,L)` computes the transferfunction of length *L*
%   of the filter defined by *g*.
%
%   See also: pfilt

complainif_argnonotinrange(nargin,2,2,mfilename);

[g,info] = comp_fourierwindow(g,L,upper(mfilename));

H=comp_transferfunction(g,L);
