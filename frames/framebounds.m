function [AF,BF]=framebounds(F);
%FRAMEBOUNDS  Frame bounds
%   Usage: fcond=framebounds(F);
%          [A,B]=framebounds(F);
%
%   `framebounds(F)` calculates the ratio $B/A$ of the frame bounds
%   of the frame given by *F*.
%
%   `[A,B]=framebounds(F)` returns the frame bounds *A* and *B* instead of
%   just their ratio.
%
%   See also: newframe, framered
  
switch(F.type)
 case {'dgt','dgtreal'}
  [AF,BF]=gabframebounds(F.g,F.a,F.M); 
 case {'dwilt','wmdct'}
  [AF,BF]=wilbounds(F.g,F.M); 
end;

if nargout<2
  % Avoid the potential warning about division by zero.
  if AF==0
    AF=Inf;
  else
    AF=BF/AF;
  end;
end;


  