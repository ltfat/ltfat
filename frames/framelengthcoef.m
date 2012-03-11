function L=framelengthcoef(coef,F);
%FRAMELENGTHCOEF  Frame length from signal
%   Usage: L=framelengthcoef(Ls,F);
%
%   `framelengthcoef(coef,F)` returns the length of the frame *F*, such that
%   *F* is long enough to expand the coefficients *coef*.
%
%   If instead a signal is given, call |framelengthsignal|_.
%
%   See also: newframe, framelengthcoef
  
switch(F.type)
 case 'dgt'
  [MN,W]=size(coef);
  L=MN/F.M*F.a;
 case 'dgtreal'
  [MN,W]=size(coef);
  L=MN/(floor(F.M/2)+1)*F.a;
 otherwise
  % handle all the bases
  L=size(coef,1);
end;
