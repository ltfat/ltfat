function outsig=frsyn(F,insig);
%FRSYN  Frame synthesis operator
%   Usage: f=frsyn(F,c);
%
%   `f=frsyn(F,c)` constructs a signal *f* from the frame coefficients *c*
%   using the frame *F*. The frame object *F* must have been created using
%   |frame|.
%
%   See also: frame, frana, plotframe
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(F)
  error('%s: First agument must be a frame definition structure.',upper(mfilename));
end;

L=framelengthcoef(F,size(insig,1));

F=frameaccel(F,L);

outsig=F.frsyn(insig);

