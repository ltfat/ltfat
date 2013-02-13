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

% Default value: works for all the bases.
red=1;

switch(F.type)
  case 'gen'
    red=size(F.g,2)/size(F.g,1);
  case 'dgt'
    red=F.M/F.a;
  case 'dgtreal'
    red=F.M/F.a;
  case {'ufilterbank','filterbank'}
    red=sum(1./F.a);
  case {'ufilterbankreal','filterbankreal'}
    red=2*sum(1./F.a);
  case 'fusion'
    red=sum(cellfun(@framered,F.frames));
 case 'tensor'
    red=prod(cellfun(@framered,F.frames));
    
end;

  