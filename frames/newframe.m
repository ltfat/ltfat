function F=newframe(ftype,varargin);
%NEWFRAME  Construct a new frame
%   Usage: F=newframe(ftype,...);
%
%   `F=newframe(ftype,...)` constructs a new frame object *F* of type
%   *ftype*. Argument following *ftype* a specific to the type of frame
%   chosen.
%
%   `newframe('dgt',g,a,M)` constructs a Gabor frame with window *g*,
%   time-shift *a* and *M* channels. See the help on |dgt|_ for more
%   information.
%
%   `newframe('dgtreal',g,a,M)` constructs a Gabor frame for real-values
%   signals with window *g*, time-shift *a* and *M* channels. See the help
%   on |dgtreal|_ for more information.
%
%   `newframe('dwilt',g,M)` constructs a Wilson basis with window *g*, and
%   *M* channels. See the help on |dwilt|_ for more information.
%
%   `newframe('wmdct',g,M)` constructs a windowed MDCT basis with window
%   *g*, and *M* channels. See the help on |wmdct|_ for more information.
%
%   Example:
%   --------
%
%   The following example create a Gabor frame for real-valued signals,
%   analysis and input signal and plots the frame coefficients:::
%
%      F=newframe('dgtreal','gauss',20,294);
%      c=framet(greasy,F);
%      plotframe(c,F);
%
%   See also: framet, iframet, plotframe
  
ftype=lower(ftype);
switch(ftype)
 case 'dgt'
  F.g=varargin{1};
  F.a=varargin{2};
  F.M=varargin{3};
 case 'dgtreal'
  F.g=varargin{1};
  F.a=varargin{2};
  F.M=varargin{3};
 case {'dwilt','wmdct'}
  F.g=varargin{1};
  F.M=varargin{2};  
end;

F.type=ftype;
F.gd={'dual',F.g};