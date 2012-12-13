function coef=framecoef2tf(F,coef);
%FRAMECOEF2TF  Convert coefficients to time-frequency plane
%   Usage: cout=framecoef2tf(F,cin);
%
%   `framecoef2tf(F,coef)` converts the frame coefficients *coef* into the
%   time-frequency plane layout. The frame object *F* must have been
%   created using |frame|_.
%
%   The time-frequency plane layout is a matrix, where the first
%   dimension indexes frequency and the second dimension time. This is
%   similar to the output format from |dgt|_ and |wmdct|_.
%
%   Not all frame types support this coefficient layout. The supported frame
%   types are: `'dgt'`, `'dgtreal'`, `'dwilt'`, `'wmdct'` and
%   `'ufilterbank'`.
%
%   See also: frame, frametf2coef, framecoef2native
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(F)
  error('%s: First agument must be a frame definition structure.',upper(mfilename));
end;

switch(F.type)
 case 'dgt'
  [MN,W]=size(coef);
  N=MN/F.M;
  coef=reshape(coef,[F.M,N,W]);  
 case 'dgtreal'
  [MN,W]=size(coef);
  M2=floor(F.M/2)+1;
  N=MN/M2;
  coef=reshape(coef,[M2,N,W]);  
 case 'dwilt'
  [MN,W]=size(coef);
  N=MN/F.M;
  coef=wil2rect(reshape(coef,[2*F.M,N/2,W]));  
 case 'wmdct'
  [MN,W]=size(coef);
  N=MN/F.M;
  coef=reshape(coef,[F.M,N,W]);  
 case 'ufilterbank'
  [MN,W]=size(coef);
  M=numel(F.gs);
  N=MN/M;
  coef=permute(reshape(coef,[N,M,W]),[2,1,3]);  
 otherwise
  error('%s: TF-plane layout not supported for this transform.',upper(mfilename));
end;

