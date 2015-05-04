function L=framelengthcoef(F,Ncoef);
%FRAMELENGTHCOEF  Frame length from coefficients
%   Usage: L=framelengthcoef(F,Ncoef);
%
%   `framelengthcoef(F,Ncoef)` returns the length of the frame *F*, such that
%   *F* is long enough to expand the coefficients of length *Ncoef*.
%
%   If instead a signal is given, call |framelength|.
%
%   See also: frame, framelength

callfun = upper(mfilename);
complainif_notenoughargs(nargin,2,callfun);
complainif_notposint(Ncoef,'Ncoef',callfun);
complainif_notvalidframeobj(F,callfun);


L = F.lengthcoef(Ncoef);
% sprintf for Octave compatibility
assert(abs(L-round(L))<1e-3,...
       sprintf('%s: There is a bug. L=%d should be an integer.',...
       upper(mfilename),L));
L=round(L);
    
% Verify the computed length
if ~(L==framelength(F,L))
    error(['%s: The coefficient number given does not correspond to a valid ' ...
           'set of coefficients for this type of frame.'],upper(mfilename));
    
end;
