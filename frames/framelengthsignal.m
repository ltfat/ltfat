function L=framelengthsignal(F,Ls);
%FRAMELENGTHSIGNAL  Frame length from signal
%   Usage: L=framelengthsignal(F,Ls);
%
%   `framelengthsignal(F,Ls)` returns the length of the frame *F*, such that
%   *F* is long enough to expand a signal of length *Ls*.
%
%   If the frame length is longer than the signal length, the signal will be
%   zero-padded by |frana|_.
%
%   If instead a set of coefficients are given, call |framelengthcoef|_.
%
%   See also: newframe, framelengthcoef
  
switch(F.type)
 case {'dgt','dgtreal'}
  L = longpar('dgt',Ls,F.a,F.M)
 case {'dwilt','wmdct'}
  L = longpar('dwilt',Ls,F.M)
 case {'filterbank','ufilterbank','filterbankreal','ufilterbankreal'}
  L = filterbanklengthsignal(F.a,Ls);
end;