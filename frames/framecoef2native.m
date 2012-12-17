function coef=framecoef2native(F,coef);
%FRAMECOEF2NATIVE  Convert coefficients to native format
%   Usage: cout=framecoef2native(F,cin);
%
%   `framecoef2native(F,coef)` converts the frame coefficients *coef* into the
%   native coefficient format of the frame. The frame object *F* must have been
%   created using |frame|_.
%
%   See also: frame, framenative2coef, framecoef2tf
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(F)
  error('%s: First agument must be a frame definition structure.',upper(mfilename));
end;

[MN,W]=size(coef);

switch(F.type)
 case 'dgt'
  N=MN/F.M;
  coef=reshape(coef,[F.M,N,W]);
  
 case 'dgtreal'
  M2=floor(F.M/2)+1;
  N=MN/M2;
  coef=reshape(coef,[M2,N,W]);
  
 case 'dwilt'
  N=MN/F.M;
  coef=reshape(coef,[2*F.M,N/2,W]);
  
 case 'wmdct'
  N=MN/F.M;
  coef=reshape(coef,[F.M,N,W]);
  
 case {'ufilterbank','ufilterbankreal'}
  N=MN/F.M;
  coef=reshape(coef,[N,F.M,W]);
  
 case {'filterbank','filterbankreal'}
  L=framelengthcoef(F,MN);
  N=L./F.a
  coef=mat2cell(coef,N,W);
    
 otherwise
  % No conversion necessary, formats are the same.
  cout=coef;
end;

