function red=framered(F);
%FRAMERED  Redundancy of a frame
%   Usage  red=framered(F);
%
%   `framered(F)` computes the redundancy of a given frame *F*. If the
%   redundancy is larger than 1 (one), the frame transform will produce more
%   coefficients than it consumes. If the redundancy is exactly 1 (one),
%   the frame is a basis.
%
%   See also: newframe, framet, framebounds
  
switch(F.type)
 case 'dgt'
  red=F.M/F.a;
 case 'dgtreal'
  red=F.M/F.a;
 case {'dwilt','wmdct'}
  red=1;
end;

  