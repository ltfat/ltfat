function outsig=frsyn(F,insig);
%FRSYN  Frame synthesis operator
%   Usage: f=frsyn(F,c);
%
%   `f=frsyn(F,c)` constructs a signal *f* from the frame coefficients *c*
%   using the frame *F*. The frame object *F* must have been created using
%   |frame|.
%
%   Examples:
%   ---------
%
%   In the following example a signal *f* is constructed through the frame
%   synthesis operator using a Gabor frame. The coefficients associated with 
%   this Gabor expansion are contained in an identity matrix. The identity 
%   matrix corresponds to a diagonal in the time-frequency plane, that is, 
%   one atom at each time position with increasing frequency.:::
%
%      a = 10;
%      M = 40;
%
%      F = frame('dgt', 'gauss', a, M);
%
%      c = framenative2coef(F, eye(40));
%
%      f = frsyn(F, c);
%
%   See also: frame, frana, plotframe
  
complainif_notenoughargs(nargin,2,'FRSYN');
complainif_notvalidframeobj(F,'FRSYN');

L=framelengthcoef(F,size(insig,1));

F=frameaccel(F,L);

outsig=F.frsyn(insig);

