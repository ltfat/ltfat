function coef=framecoef2native(F,coef);
%FRAMECOEF2NATIVE  Convert coefficients to native format
%   Usage: cout=framecoef2native(F,cin);
%
%   `framecoef2native(F,coef)` converts the frame coefficients *coef* into the
%   native coefficient format of the frame. The frame object *F* must have been
%   created using |frame|.
%
%   See also: frame, framenative2coef, framecoef2tf
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(F)
  error('%s: First agument must be a frame definition structure.',upper(mfilename));
end;

[MN,W]=size(coef);

if isfield(F,'coef2native')
    coef=F.coef2native(coef,size(coef));
else
    
    switch(F.type)                                
      case {'filterbank','filterbankreal'}
        L=framelengthcoef(F,MN);
        N=L./F.a
        coef=mat2cell(coef,N,W);
  
    end;
end;
