function [g,L,info] = nonsepgabpars_from_window(g,a,M,s,L,callfun)
%NONSEPGABPARS_FROM_WINDOW  Compute g and L from window
%   Usage: [g,g.info,L] = gabpars_from_window(f,g,a,M);
%
%   Use this function if you know a window and a lattice
%   for the NONSEPDGT. The function will calculate a transform length L and
%   evaluate the window g into numerical form.
%
%   If the transform length is unknown (as it usually is unless explicitly
%   specified by the user), set L to be [] in the input to this function.
  
if nargin<6
  stacknames=dbstack;  
  callfun=stacknames(2).name;
end;

if isempty(L)
  if isnumeric(g)
    L=length(g);
  else
    L=nonsepdgtlengthsignal(1,a,M,s);
  end;
else
  Lcheck=nonsepdgtlengthsignal(L,a,M,s);
  if Lcheck~=L
    error('%s: Invalid transform size L',upper(mfilename));
  end;
end;

[g,info] = comp_window(g,a,M,L,s,'NONSEPGABDUAL');

if (info.isfir)  
  if info.istight
    g=g/sqrt(2);
  end;  
end;
