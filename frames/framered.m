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
%     F=newframe('dgt','gauss','dual',30,40);
%     framered(F)
%
%   The redundancy of a basis is always one:::
%
%     F=newframe('wmdct','gauss','dual',40);
%     framered(F)
%
%   See also: newframe, frana, framebounds

% Default value: works for all the bases.
red=1;

switch(F.type)
 case 'gen'
  red=size(F.ga,2)/size(F.ga,1);
 case 'dgt'
  red=F.M/F.a;
 case 'dgtreal'
  red=F.M/F.a;
end;

  