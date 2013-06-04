function coef=framenative2coef(F,coef);
%FRAMENATIVE2COEF  Convert coefficient from native format
%   Usage: cout=framenative2coef(F,cin);
%
%   `framenative2coef(F,coef)` convert frame coefficients from the native
%   format of the transform into the common column format.
%
%   See also: frame, framecoef2native
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(F)
  error('%s: First agument must be a frame definition structure.',upper(mfilename));
end;

coef=F.native2coef(coef);
