function [coef]=comp_dwilt(f,g,M)
%COMP_DWILT  Compute Discrete Wilson transform.
%   

L=size(f,1);
Lwindow=size(g,1);

if Lwindow<L
  coef=comp_dwilt_fb(f,g,M,L);
else
  coef=comp_dwilt_long(f,g,M,L);
end;


