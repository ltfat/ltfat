function ftype=frametype(F);
%FRAMETYPE  Frame type
%   Usage: L=frametype(F);
%
%   `frametype(F)` returns a text string describing the type of the
%   frame. The string is the same as the one used to create the frame in
%   |newframe|_.
%
%   See also: newframe
  
ftype=F.type;