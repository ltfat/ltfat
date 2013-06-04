function L=framelength(F,Ls);
%FRAMELENGTH  Frame length from signal
%   Usage: L=framelength(F,Ls);
%
%   `framelength(F,Ls)` returns the length of the frame *F*, such that
%   *F* is long enough to expand a signal of length *Ls*.
%
%   If the frame length is longer than the signal length, the signal will be
%   zero-padded by |frana|.
%
%   If instead a set of coefficients are given, call |framelengthcoef|.
%
%   See also: frame, framelengthcoef
  
L=F.length(Ls);
