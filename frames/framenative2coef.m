function coef=framenative2coef(F,coef);
%FRAMENATIVE2COEF  Convert coefficient from native format
%   Usage: coef=framenative2coef(F,coef);
%
%   `framenative2coef(F,coef)` converts the frame coefficients from the 
%   native format of the transform into the common column format.
%
%   See also: frame, framecoef2native
  
complainif_notenoughargs(nargin,2,'FRAMENATIVE2COEF');
complainif_notvalidframeobj(F,'FRAMENATIVE2COEF');

% .native2coef is not a mandatory field
if isfield(F,'native2coef')
   coef=F.native2coef(coef);
end
