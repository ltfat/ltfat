function coef=frametf2coef(F,coef)
%FRAMETF2COEF  Convert coefficients from TF-plane format
%   Usage: coef=frametf2coef(F,coef);
%
%   `frametf2coef(F,coef)` converts the frame coefficients from the
%   time-frequency plane layout into the common column format.
%
%   See also: frame, framecoef2tf, framecoef2native

complainif_notenoughargs(nargin,2,'FRAMETF2COEF');
complainif_notvalidframeobj(F,'FRAMETF2COEF');


switch(F.type)
 case {'dgt','dgtreal','wmdct'}
  [M,N,W]=size(coef);
  coef=reshape(coef,[M*N,W]);
 case {'dwilt'}
  coef=framenative2coef(F,rect2wil(coef));
 case {'ufilterbank'}
   coef=permute(coef,[2,1,3]);
   [M,N,W]=size(coef);
   coef=reshape(coef,[M*N,W]);
 case {'ufwt','uwfbt','uwpfbt'}
  coef = F.native2coef(permute(coef,[2,1,3])); 
 otherwise
  error('%s: TF-plane layout not supported for this transform.',upper(mfilename));
end;



