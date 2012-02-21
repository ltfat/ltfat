function plotdgt(coef,a,varargin)
%PLOTDGT  Plot DGT coefficients
%   Usage: plotdgt(coef,a);
%          plotdgt(coef,a,fs);
%          plotdgt(coef,a,fs,dynrange);
%
%   `plotdgt(coef,a)` plots the Gabor coefficients *coef*. The coefficients
%   must have been produced with a timeshift of *a*.
%
%   `plotdgt(coef,a,fs)` does the same assuming a sampling rate of
%   *fs* Hz of the original signal.
%
%   `plotdgt(coef,a,fs,dynrange)` additionally limits the dynamic range.
%   
%   `plotdgt` supports all the optional parameters of |tfplot|_. Please see
%   the help of |tfplot|_ for an exhaustive list.
%
%   See also:  dgt, tfplot, sgram, plotdgtreal

%   AUTHOR : Peter Soendergaard.
%   TESTING: NA
%   REFERENCE: NA

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'ltfattranslate','tfplot'};

[flags,kv,fs]=ltfatarghelper({'fs','dynrange'},definput,varargin);

M=size(coef,1);

% Move zero frequency to the center and Nyquest frequency to the top.
if rem(M,2)==0
  coef=circshift(coef,M/2-1);
  yr=[-1+2/M, 1];
else
  coef=circshift(coef,(M-1)/2);
  yr=[-1+2/M, 1-2/M];
end;

tfplot(coef,a,yr,'argimport',flags,kv);

