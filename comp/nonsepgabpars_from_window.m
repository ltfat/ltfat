function [g,L,info] = nonsepgabpars_from_window(g,a,M,lt,L,callfun)
%NONSEPGABPARS_FROM_WINDOW  Compute g and L from window
%   Usage: [g,g.info,L] = gabpars_from_window(f,g,a,M,lt,L);
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
    L=dgtlength(1,a,M,lt);
  end;
else
  Lcheck=dgtlength(L,a,M,lt);
  if Lcheck~=L
    error('%s: Invalid transform size L',upper(mfilename));
  end;
end;

[g,info] = comp_window(g,a,M,L,lt,'NONSEPGABDUAL');

if (info.isfir)  
  if info.istight
    g=g/sqrt(2);
  end;  
end;
