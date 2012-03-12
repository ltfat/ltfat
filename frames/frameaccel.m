function F=frameaccel(F,L);  
%FRAMEACCEL  Precompute structures
%   Usage: F=frameaccel(F,L);
%
%   `F=frameaccel(F,L)` precomputes certain structures that makes the basic
%   frame operations |framet|_ and |iframet|_ faster (like instantiating the
%   window from a textual description). If you only need to call the
%   routines once, calling `frameaccel` first will not provide any total
%   gain, but if you are repeatedly calling these routines, for instance in
%   an iterative algorithm, is will be a benefit.
%
%   Notice that you need to input the frame length *L*, so this routines
%   is only a benefit if *L* stays fixed.
%
%   See also: newframe, framet, framelengthsignal, framelengthcoef
  
switch(F.type)
 case {'dgt','dgtreal'}
  [F.g, F.g_info]  = gabwin(F.g,F.a,F.M,L);
  [F.gd,F.gd_info] = gabwin(F.gd,F.a,F.M,L);
 case {'dwilt','wmdct'}
  [F.g, F.g_info]  = wilwin(F.g,F.M,L);
  [F.gd,F.gd_info] = wilwin(F.gd,F.M,L);
end;

  