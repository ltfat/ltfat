function plotdgt(coef,a,varargin)
%DGTPLOT  Plot DGT coefficients.
%   Usage: plotdgt(coef,a);
%          plotdgt(coef,a,fs);
%          plotdgt(coef,a,fs,dynrange);
%
%   PLOTDGT(coef,a) will plot the Gabor coefficient coefficients
%   coef. The coefficient must have been produce with a timeshift of _a.
%
%   PLOTDGT(coef,a,fs) will do the same assuming a sampling rate of
%   fs Hz of the original signal.
%
%   PLOTDGT(coef,a,fs,dynrange) will additionally limit the dynamic
%   range.
%   
%   PLOTDGT supports all the optional parameters of TFPLOT. Please
%   see the help of TFPLOT for an exhaustive list.
%
%   See also:  dgt, tfplot, sgram, plotdgtreal

%   AUTHOR : Peter Soendergaard.
%   TESTING: NA
%   REFERENCE: NA

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'tfplot'};

[flags,kv,fs]=ltfatarghelper({'fs','dynrange'},definput,varargin);

M=size(coef,1);

% Move zero frequency to the center and Nyquest frequency to the top.
if rem(M,2)==0
  coef=circshift(coef,M/2-1);
  yr=(-M/2+1:M/2)/(M/2);
else
  coef=circshift(coef,(M-1)/2);
  yr=(-(M-1)/2:(M-1)/2)/((M-1)/2);
end;

if ~isempty(kv.fs)
  yr=yr*fs/2;
end;


tfplot(coef,a,yr,'argimport',flags,kv);