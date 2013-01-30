function c = fwt(f,h,J,varargin)
%FWT   Fast Wavelet Transform 
%   Usage:  c = fwt(f,h,J);
%           c = fwt(f,h,J,...);
%
%   Input parameters:
%         f     : Input data.
%         h     : Analysis Wavelet Filterbank. 
%         J     : Number of filterbank iterations.
%         a     : Explicitly defined subsampling factors.
%
%   Output parameters:
%         c      : Coefficients stored in a cell-array.
%
%   `c=fwt(f,h,J)` returns wavelet coefficients *c* of the input signal *f*
%   using *J* iterations of the basic wavelet filterbank defined by *h*.
%   The fast wavelet transform algorithm (or Mallat's algorithm) is employed.
%   If *f* is a matrix, the transformation is applied to each of *W* columns.
%   
%   The coefficents *c* are Discrete Wavelet transform of the input signal *f*,
%   if *h* defines two-channel wavelet filterbank. The function can apply
%   the Mallat's algorithm using basic filter banks with any number of the
%   channels. In such case, the transform have different name.    
%
%   The basic analysis wavelet filterbank $h$ can be passed in several formats. 
%   The formats are the same as for the |fwtinit|_ function. The simplest
%   is passing a cell array, whose first element is the name of the function
%   defining the basic wavelet filters (|wfilt_|_ prefix) and the other elements are 
%   the parameters of the function. e.g. `{'db',10}` calls  `wfilt_db(10)` internally.
%   The second possible format of $h$ is to pass cell array of one dimensional
%   numerical vectors directly defining the wavelet filter impulse responses.
%   In this case, outputs of the filters are subsampled by a factor equal
%   to the number of the filters. For creating completely custom filterbanks use the
%   |fwtinit|_ function.
%   The third option is to pass structure obtained from the |fwtinit|_
%   function.
%   
%   The coefficients in the cell array $c\{jj\}$ for $jj=1,\ldots,J\cdot(N-1)+1$,
%   where $N$ is number of filters in the basic wavelet filterbank, are ordered
%   with inceasing central frequency of the corresponding effective filter frequency
%   responses or equivalently with decreasing wavelet scale.
%
%   If the input *f* is matrix with *W* columns, each element of the cell
%   array $c\{jj\}$ is a matrix with *W* columns with coefficients belonging
%   to the appropriate input channel.
%
%   The following flag groups are supported:
%
%         'per','zpd','sym','symw','asym','asymw','ppd','sp0'
%                Type of the boundary handling.
%
%
%   Boundary handling:
%   ------------------
%
%   The default periodic extension considers the input signal as it was a one
%   period of some infinite periodic signal as is natural for the transforms
%   based on the FFT.  The resulting wavelet representation is
%   non-expansive, that is if the input signal length is a multiple of a
%   $J$-th power of the subsampling factor, the total number of coefficients
%   is equal to the input signal length. 
%
%   If the input signal length is not a multiple of a $J$-th power of the
%   subsampling factor, the processed signal is padded internally by
%   repeating the last sample at each step of the transform to the next
%   multiple of the subsampling factor rather than doing the prior explicit
%   padding.  In addition, the periodic extension restrict the input signal
%   length to be greater than a certain length.
%
%   `fwt(f,h,J,ext)` with `ext` other than `'per'` computes a slightly
%   redundant wavelet representation of the input signal *f* with the chosen
%   boundary extension *ext*.
%
%   The default periodic extension can result in "false" high wavelet
%   coefficients near the boundaries due to the possible discontinuity
%   introduced by the periodic extension. The custom extension can diminish
%   this phenomenon. The extensions are done at each level of the transform
%   internally rather than doing the prior explicit padding. The supported
%   possibilities are:
%
%     * `'zpd'` - zeros are considered outside of the signal (coefficient) support. 
%
%     * `'sym'` - half-point symmetric extension.
%
%     * `'symw'` - whole-point symmetric extension
%
%     * `'asym'` - half-point antisymmetric extension
%
%     * `'asymw'` - whole point antisymmetric extension
%
%     * `'ppd'` - periodic padding, same as `'ppd'` but the result is expansive representation
%
%     * `'sp0'` - repeating boundary sample
%
%   Note that the same flag have to be used in the call of the inverse transform
%   function `ifwt`.
%
%   Examples:
%   ---------
%   
%   A simple example of calling the |fwt|_ function:::
% 
%     f = gspi;
%     J = 10;
%     c = fwt(f,{'db',8},J);
%     plotfwt(c);
%
%   See also: ifwt, fwtinit
%
%   References: ma98  


if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(J) || ~isscalar(J)
  error('%s: "J" must be a scalar.',upper(mfilename));
end;

if(J<1 || rem(J,1)~=0)
   error('%s: J must be a positive integer.',upper(mfilename)); 
end

% Initialize the wavelet filters structure
h = fwtinit(h,'ana');

%% ----- step 0 : Check inputs -------
definput.import = {'fwt'};
[flags,kv]=ltfatarghelper({},definput,varargin);


%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
% TO DO: if elements of h.a are not equal and the flag is 'per',
% do zero padding to the next multiple of 2^J
[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,upper(mfilename),0);
if(Ls<2)
   error('%s: Input signal seems not to be a vector of length > 1.',upper(mfilename));  
end


%% ----- step 2 : Check whether the input signal is long enough
% input signal length is not restricted for expansive wavelet transform (extension type other than the default 'per')
flen = length(h.filts{1}.h);
if(strcmp(flags.ext,'per'))
     minLs = (h.a(1)^J-1)*(flen-1); % length of the longest equivalent filter -1
   if Ls<minLs
     error('%s: Input signal length is %d. Minimum signal length is %d or use %s flag instead. \n',upper(mfilename),Ls,minLs,'''ppd''');
   end;
end


%% ----- step 3 : Run computation
 c = comp_fwt_all(f,h.filts,J,h.a,'dec',flags.ext);


