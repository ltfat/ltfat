function [AF,BF]=framebounds(F,varargin);
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
%   `framebounds(F,'s')` returns the framebounds of the synthesis frame
%   instead of those of the analysis frame.
%
%   See also: newframe, framered

definput.keyvals.L=[];
definput.flags.system={'a','s'};
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);
 
% Default values, works for the pure frequency transforms.
AF=1;
BF=1;

if isfield(F,'ga')
  if flags.do_a
    g=F.ga;
  else
    g=F.gs;
  end;
end;

switch(F.type)
 case 'gen'
  V=svd(g);
  AF=min(V)^2;
  BF=max(V)^2;
 case {'dgt','dgtreal'}
  [AF,BF]=gabframebounds(g,F.a,F.M); 
 case {'dwilt','wmdct'}
  [AF,BF]=wilbounds(g,F.M); 
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


  