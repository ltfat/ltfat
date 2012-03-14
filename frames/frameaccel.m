function F=frameaccel(F,L);  
%FRAMEACCEL  Precompute structures
%   Usage: F=frameaccel(F,L);
%
%   `F=frameaccel(F,L)` precomputes certain structures that makes the basic
%   frame operations |frana|_ and |frsyn|_ faster (like instantiating the
%   window from a textual description). If you only need to call the
%   routines once, calling `frameaccel` first will not provide any total
%   gain, but if you are repeatedly calling these routines, for instance in
%   an iterative algorithm, is will be a benefit.
%
%   Notice that you need to input the frame length *L*, so this routines
%   is only a benefit if *L* stays fixed.
%
%   See also: newframe, frana, framelengthsignal, framelengthcoef
  
switch(F.type)
 case {'dgt','dgtreal'}
  [F.ga,F.g_info]  = gabwin(F.ga,F.a,F.M,L);
  [F.gs,F.gs_info] = gabwin(F.gs,F.a,F.M,L);
 case {'dwilt','wmdct'}
  [F.ga,F.g_info]  = wilwin(F.ga,F.M,L);
  [F.gs,F.gs_info] = wilwin(F.gs,F.M,L);
end;

  