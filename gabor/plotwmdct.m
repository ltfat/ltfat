function plotwmdct(coef,varargin)
%PLOTWMDCT  Plot WMDCT coefficients.
%   Usage: plotwmdct(coef,);
%          plotwmdct(coef,fs);
%          plotwmdct(coef,fs,dynrange);
%
%   PLOTWMDCT(coef) will plot coefficients from WMDCT.
%
%   PLOTWMDCT(coef,fs) will do the same assuming a sampling rate of
%   fs Hz of the original signal.
%
%   PLOTWMDCT(coef,fs,dynrange) will additionally limit the dynamic
%   range.
%   
%   PLOTWMDCT supports all the optional parameters of TFPLOT. Please
%   see the help of TFPLOT for an exhaustive list.
%
%   See also:  wmdct, tfplot, sgram, plotdgt

%   AUTHOR : Peter Soendergaard.
%   TESTING: NA
%   REFERENCE: NA

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'ltfattranslate','tfplot'};

[flags,kv,fs]=ltfatarghelper({'fs','dynrange'},definput,varargin);

M=size(coef,1);

yr=[.5/M, 1-.5/M];

tfplot(coef,M,yr,'argimport',flags,kv);