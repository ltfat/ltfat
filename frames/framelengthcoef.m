function L=framelengthcoef(F,Ncoef);
%FRAMELENGTHCOEF  Frame length from coefficients
%   Usage: L=framelengthcoef(F,Ncoef);
%
%   `framelengthcoef(F,Ncoef)` returns the length of the frame *F*, such that
%   *F* is long enough to expand the coefficients of length *Ncoef*.
%
%   If instead a signal is given, call |framelengthsignal|_.
%
%   See also: newframe, framelengthsignal
  
switch(F.type)
 case 'dgt'
  L=Ncoef/F.M*F.a;
 case 'dgtreal'
  L=Ncoef/(floor(F.M/2)+1)*F.a;
 otherwise
  % handle all the bases
  L=Ncoef;
end;
