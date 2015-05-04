function red=framered(F);
%FRAMERED  Redundancy of a frame
%   Usage  red=framered(F);
%
%   `framered(F)` computes the redundancy of a given frame *F*. If the
%   redundancy is larger than 1 (one), the frame transform will produce more
%   coefficients than it consumes. If the redundancy is exactly 1 (one),
%   the frame is a basis.
%
%   Examples:
%   ---------
%
%   The following simple example shows how to obtain the redundancy of a
%   Gabor frame:::
%
%     F=frame('dgt','gauss',30,40);
%     framered(F)
%
%   The redundancy of a basis is always one:::
%
%     F=frame('wmdct','gauss',40);
%     framered(F)
%
%   See also: frame, frana, framebounds

complainif_notenoughargs(nargin,1,'FRAMERED');
complainif_notvalidframeobj(F,'FRAMERED');

% .red field is mandatory so no checking here
red=F.red;
  
