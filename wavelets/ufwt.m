function c = ufwt(f,h,J,varargin)
%UFWT  Undecimated Fast Wavelet Transform 
%   Usage:  c = ufwt(f,h,J);
%           c = ufwt(f,h,J,...);
%
%   Input parameters:
%         f     : Input data.
%         h     : Analysis Wavelet Filterbank. 
%         J     : Number of filterbank iterations.
%
%   Output parameters:
%         c      : Coefficients stored in a cell-array.
%
%   `c=ufwt(f,h,J)` computes redundant time (or shift) invariant
%   wavelet representation *c* of the input signal *f* using the "a-trous"
%   algorithm.

%   The coefficents *c* are so coalled Undecimated Discrete Wavelet transform
%   of the input signal *f*, if *h* defines two-channel wavelet filterbank.
%   Other names for this version of the wavelet transform are: the
%   time-invariant wavelet transform, the stationary wavelet transform, 
%   maximal overlap discrete wavelet transform or even the "continuous" 
%   wavelet transform (as the time step is one sample).
%   However, the function accepts any number channels (further referred to as $N$)
%   in the basic wavelet filterbank. 
%   By default, the redundancy of the representanion is $J*(N-1)+1$ and
%   `length(c{jj})=length(f)` for all *jj*.
%
%   For all accepted formats of the parameter *h* see the |fwt|_ function.
%
%   The coefficients in the cell array $c\{j\}$ for $j=1,\ldots,J\cdot(N-1)+1$ are ordered
%   with inceasing central frequency of the corresponding effective filter frequency
%   response or equivalently with decreasing wavelet scale.
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
%   `c=ufwt(f,h,J,ext)` with `ext` other than `'per'` computes even more
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
%   Examples:
%   ---------
%   
%   A simple example of calling the |ufwt|_ function:::
% 
%     f = gspi;
%     J = 10;
%     c = ufwt(f,{'db',8},J);
%     plotfwt(c);
%
%   See also: iunfwt, fwtinit
%
%   References: ma98  



if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(J) || ~isscalar(J)
  error('%s: "J" must be a scalar.',upper(mfilename));
end;

if(J<1 && rem(a,1)~=0)
   error('%s: J must be a positive integer.',upper(mfilename)); 
end

% Initialize the wavelet filters structure
h = fwtinit(h);

%% ----- step 0 : Check inputs -------
definput.import = {'fwt'};
[flags,kv]=ltfatarghelper({},definput,varargin);


%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,upper(mfilename),0);
if(Ls<2)
   error('%s: Input signal seems not to be a vector of length > 1.',upper(mfilename));  
end
 


%% ----- step 2 : Check whether the input signal is long enough
% input signal length is not restricted for expansive wavelet transform (extension type other than the default 'per')
flen = length(h{1});
if(strcmp(flags.ext,'per'))
     minLs = (a(1)^(J-1))*(flen-1); % length of the longest upsampled filter - 1
   if Ls<minLs
     error('%s: Input signal length is %d. Minimum signal length is %d or use %s flag instead. \n',upper(mfilename),Ls,minLs,'''ppd''');
   end;
end



%% ----- step 3 : Run computation
 c = comp_fwt_all(f,h,J,a,'undec',flags.ext);


