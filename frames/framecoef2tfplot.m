function coef=framecoef2tfplot(F,coef)
%FRAMECOEF2TFPLOT  Convert coefficients to time-frequency plane matrix
%   Usage: cout=framecoef2tfplot(F,cin);
%
%   `framecoef2tfplot(F,coef)` converts the frame coefficients *coef* into
%   the time-frequency plane layout matrix. The frame object *F* must have 
%   been created using |frame|. The function acts exactly as 
%   |framecoef2tf| for frames which admit regular (rectangular) sampling
%   of a time-frequency plane and converts irregularly sampled coefficients
%   to a rectangular matrix. This is usefull for custom plotting.
%
%   See also: frame, frametf2coef, framecoef2native, blockplot
  
complainif_notenoughargs(nargin,2,'FRAMECOEF2TFPLOT');
complainif_notvalidframeobj(F,'FRAMECOEF2TFPLOT');

switch(F.type)
   % We could have done a try-catch block here, but it is slow
   case {'dgt','dgtreal','dwilt','wmdct','ufilterbank','ufwt','uwfbt','uwpfbt'} 
    coef=framecoef2tf(F,coef);
    return;
   case 'fwt'
    coef = comp_fwtpack2cell(F,coef);
   case {'wfbt','wpfbt','filterbank','filterbankreal'}
    coef = F.coef2native(coef,size(coef));   
end

switch(F.type)
 case {'fwt','wfbt','wpfbt','filterbank','filterbankreal'}
  coef = comp_cellcoef2tf(coef);
 otherwise
  error('%s: TF-plane plot not supported for this transform.',upper(mfilename));
end;

