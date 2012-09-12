function L=framelength(F,Ls);
%FRAMELENGTH  Frame length from signal
%   Usage: L=framelength(F,Ls);
%
%   `framelength(F,Ls)` returns the length of the frame *F*, such that
%   *F* is long enough to expand a signal of length *Ls*.
%
%   If the frame length is longer than the signal length, the signal will be
%   zero-padded by |frana|_.
%
%   If instead a set of coefficients are given, call |framelengthcoef|_.
%
%   See also: newframe, framelengthcoef
  
% Default value, the frame works for all input lengths
L=Ls;
  
switch(F.type)
 case {'dgt','dgtreal'}
  L = dgtlength(Ls,F.a,F.M);
 case {'dwilt','wmdct'}
  L = longpar('dwilt',Ls,F.M);
 case {'gen'}
  L = size(F.ga,1);
 case {'filterbank','ufilterbank','filterbankreal','ufilterbankreal'}
  L = filterbanklength(Ls,F.a);
end;