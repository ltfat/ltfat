function coef=plotdgtreal(coef,a,M,varargin)
%PLOTDGTREAL  Plot DGTREAL coefficients
%   Usage: plotdgtreal(coef,a,M);
%          plotdgtreal(coef,a,M,fs);
%          plotdgtreal(coef,a,M,fs,dynrange);
%
%   `plotdgtreal(coef,a,M)` plots Gabor coefficient from |dgtreal|. The
%   parameters *a* and *M* must match those from the call to |dgtreal|.
%
%   `plotdgtreal(coef,a,M,fs)` does the same assuming a sampling rate of *fs*
%   Hz of the original signal.
%
%   `plotdgtreal(coef,a,M,fs,dynrange)` additionally limits the dynamic
%   range.
%
%   `C=plotdgtreal(...)` returns the processed image data used in the
%   plotting. Inputting this data directly to `imagesc` or similar
%   functions will create the plot. This is usefull for custom
%   post-processing of the image data.
%
%   `plotdgtreal` supports all the optional parameters of |tfplot|. Please
%   see the help of |tfplot| for an exhaustive list.
%
%   See also:  dgtreal, tfplot, sgram, plotdgt

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: NA
%   REFERENCE: NA

complainif_notenoughargs(nargin,3,mfilename);
complainif_notposint(a,'a',mfilename);
complainif_notposint(M,'M',mfilename);

definput.import={'ltfattranslate','tfplot'};

[flags,kv,fs]=ltfatarghelper({'fs','dynrange'},definput,varargin);

if rem(M,2)==0
  yr=[0,1];
else
  yr=[0,1-2/M];
end;

coef=tfplot(coef,a,yr,'argimport',flags,kv);

if nargout<1
    clear coef;
end
