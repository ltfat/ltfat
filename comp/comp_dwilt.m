function [coef]=comp_dwilt(f,g,M,L)
%COMP_DWILT  Compute Discrete Wilson transform.
%   

% XXX This function should not take L as an input

Lwindow=size(g,1);

if Lwindow<L
  coef=comp_dwilt_fb(f,g,M,L);
else
  coef=comp_dwilt_long(f,g,M,L);
end;


