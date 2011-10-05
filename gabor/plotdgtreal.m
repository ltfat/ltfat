function plotdgtreal(coef,a,M,varargin)
%PLOTDGTREAL  Plot DGTREAL coefficients.
%   Usage: plotdgtreal(coef,a,M);
%          plotdgtreal(coef,a,M,fs);
%          plotdgtreal(coef,a,M,fs,dynrange);
%
%   PLOTDGTREAL(coef,a,M) will plot Gabor coefficient from DGTREAL. The
%   parameters _a and M must match those from the call to DGTREAL.
%
%   PLOTDGTREAL(coef,a,M,fs) will do the same assuming a sampling rate of
%   fs Hz of the original signal.
%
%   PLOTDGTREAL(coef,a,M,fs,dynrange) will additionally limit the dynamic
%   range.
%   
%   PLOTDGTREAL supports all the optional parameters of TFPLOT. Please
%   see the help of TFPLOT for an exhaustive list.
%
%   See also:  dgtreal, tfplot, sgram, plotdgt

%   AUTHOR : Peter Soendergaard.
%   TESTING: NA
%   REFERENCE: NA

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'ltfattranslate','tfplot'};

[flags,kv,fs]=ltfatarghelper({'fs','dynrange'},definput,varargin);

if rem(M,2)==0
  yr=[0,1];
else
  yr=[0,1-2/M];
end;

tfplot(coef,a,yr,'argimport',flags,kv);