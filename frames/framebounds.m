function [AF,BF]=framebounds(F,L);
%FRAMEBOUNDS  Frame bounds
%   Usage: fcond=framebounds(F);
%          [A,B]=framebounds(F);
%
%   `framebounds(F)` calculates the ratio $B/A$ of the frame bounds
%   of the frame given by *F*.
%
%   `framebounds(F,L)` additionally specifies the length of the
%   transform. This is necessary for the `'fft'` frame.
%
%   `[A,B]=framebounds(F)` returns the frame bounds *A* and *B* instead of
%   just their ratio.
%
%   See also: newframe, framered

% Default values, works for the pure frequency transforms.
AF=1;
BF=1;

switch(F.type)
 case 'gen'
  V=svd(F.ga);
  AF=min(V)^2;
  BF=max(V)^2;
 case {'dgt','dgtreal'}
  [AF,BF]=gabframebounds(F.ga,F.a,F.M); 
 case {'dwilt','wmdct'}
  [AF,BF]=wilbounds(F.ga,F.M); 
 case 'fft'
  AF=L;
  BF=L;
 case 'fftreal'
  AF=F.L;
  BF=F.L;  
end;

if nargout<2
  % Avoid the potential warning about division by zero.
  if AF==0
    AF=Inf;
  else
    AF=BF/AF;
  end;
end;


  