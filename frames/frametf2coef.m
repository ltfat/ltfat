function coef=frametf2coef(F,coef);
%FRAMETF2COEF  Convert coefficients from TF-plane format
%   Usage: cout=frametf2coef(F,cin);
%
%   `frametf2coef(F,coef)` convert frame coefficients from the
%   time-frequency plane layout into the common column format.
%
%   See also: frame, framecoef2tf, framecoef2native
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(F)
  error('%s: First agument must be a frame definition structure.',upper(mfilename));
end;

switch(F.type)
 case {'dgt','dgtreal','wmdct'}
  [M,N,W]=size(coef);
  coef=reshape(coef,[M*N,W]); 
 case {'wmdct','ufilterbank'}
  coef=framenative2coef(F,rect2wil(coef));
 case {'ufilterbank'}
  coef=permute(coef,[2,1,3]);
  [M,N,W]=size(coef);
  coef=reshape(coef,[M*N,W]); 
 otherwise
  error('%s: TF-plane layout not supported for this transform.',upper(mfilename));
end;

