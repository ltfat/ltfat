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
%   See also: frame, framelengthcoef, frameclength
  
callfun = upper(mfilename);
complainif_notenoughargs(nargin,2,callfun);
complainif_notposint(Ls,'Ls',callfun);
complainif_notvalidframeobj(F,callfun);

% .length field is mandatory
L=F.length(Ls);
