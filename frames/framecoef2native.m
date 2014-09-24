function coef=framecoef2native(F,coef)
%FRAMECOEF2NATIVE  Convert coefficients to native format
%   Usage: coef=framecoef2native(F,coef);
%
%   `framecoef2native(F,coef)` converts the frame coefficients *coef* into 
%   the native coefficient format of the frame. The frame object *F* must 
%   have been created using |frame|.
%
%   See also: frame, framenative2coef, framecoef2tf
  
complainif_notenoughargs(nargin,2,'FRAMECOEF2NATIVE');
complainif_notvalidframeobj(F,'FRAMECOEF2NATIVE');

[MN,W]=size(coef);

% .coef2native field is not mandatory since for some frames, both
% coefficient formats are identical
if isfield(F,'coef2native')
    coef=F.coef2native(coef,size(coef));
end;
