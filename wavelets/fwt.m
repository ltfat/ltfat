function [c,Lc] = fwt(f,h,J,varargin)
%FWT   Fast Wavelet Transform 
%   Usage:  c = fwt(f,h,J);
%           c = fwt(f,h,J,'ext',ext,'dim',dim);
%
%   Input parameters:
%         f     : Input data.
%         h     : Analysis Wavelet Filterbank. 
%         J     : Number of filterbank iterations.
%
%   Output parameters:
%         c      : Coefficients stored as matrix.
%
%   `c=fwt(f,h,J)` returns wavelet coefficients *c* of the input signal *f*
%   using *J* iterations of the basic wavelet filterbank defined by *h*.
%   The fast wavelet transform algorithm (Mallat's algorithm) is employed.
%   
%   The coefficents *c* are Discrete Wavelet transform of the input signal *f*,
%   if *h* defines two-channel wavelet filterbank. The function can apply
%   the Mallat's algorithm using basic filter banks with any number of the
%   channels. In such case, the transform have a different name.    
%
%   The basic analysis wavelet filterbank $h$ can be passed in several formats. 
%   The formats are the same as for the |fwtinit|_ function.
%
%     1) The simplest is passing a cell array, whose first element is the
%        name of the function defining the basic wavelet filters (`wfilt_`
%        prefix) and the other elements are the parameters of the
%        function. e.g. `{'db',10}` calls `wfilt_db(10)` internally.
%
%     2) The second possible format of $h$ is to pass cell array of one
%        dimensional numerical vectors directly defining the wavelet filter
%        impulse responses.  In this case, outputs of the filters are
%        subsampled by a factor equal to the number of the filters. For
%        creating completely custom filterbanks use the |fwtinit|_ function.
%
%     3) The third option is to pass a structure obtained from the
%        |fwtinit|_ function.
%   
%   If *f* is row/collumn vector, the subbands *c* are stored in a single
%   collumn/row in a consecutive order with respect to the inceasing central frequency
%   of the corresponding effective filter frequency responses or
%   equivalently with decreasing wavelet scale. The lengths of subbands are
%   stored in *Lc*.
%
%   If the input *f* is matrix with *W* columns, the transform is applied to
%   each collumn (*dim=1*) or row (*dim=2*) and the dimension is preserved
%   in the output coefficients.
%
%   Boundary handling:
%   ------------------
%
%   `fwt(f,h,J,'per')` (default) uses the periodic extension which considers
%   the input signal as it was a one period of some infinite periodic signal
%   as is natural for transforms based on the FFT. The resulting wavelet
%   representation is non-expansive, that is if the input signal length is a
%   multiple of a $J$-th power of the subsampling factor and the filterbank
%   is critically subsampled, the total number of coefficients is equal to
%   the input signal length. The input signal is padded with zeros to the
%   next legal length internally.
%   
%   The default periodic extension can result in "false" high wavelet
%   coefficients near the boundaries due to the possible discontinuity
%   introduced by the zero padding and periodic boundary treatment.
%
%   `fwt(f,h,J,ext)` with `ext` other than `'per'` computes a slightly
%   redundant wavelet representation of the input signal *f* with the chosen
%   boundary extension *ext*. The redundancy (expansivity) of the
%   represenation is the price to pay for using general filterbank and
%   custom boundary treatment.  The extensions are done at each level of the
%   transform internally rather than doing the prior explicit padding.
%   
%
%   The supported possibilities are:
%
%     'per'    Zero padding to the next legal length and periodic boundary
%              extension. This is the default.
%
%     'zero'   Zeros are considered outside of the signal (coefficient)
%              support. 
%
%     'even'   Even symmetric extension.
%
%     'odd'    Odd symmetric extension.
%
%   Note that the same flag has to be used in the call of the inverse transform
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
%     %plotfwt(c,44100,90);
%
%   See also: ifwt, fwtinit, wavpack2cell, plotfwt
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
definput.keyvals.dim = [];
[flags,kv,dim]=ltfatarghelper({'dim'},definput,varargin);

%% ----- step 1 : Verify f and determine its length -------
[f,L,Ls,~,dim,~,order]=assert_sigreshape_pre(f,[],dim,upper(mfilename));

% Determine next legal input data length.
L = fwtlength(Ls,h,J,flags.ext);

% Pad with zeros if the safe length L differ from the Ls.
if(Ls~=L)
   f=postpad(f,L); 
end

%% ----- step 2 : Determine number of ceoefficients in each subband
Lc = fwtclength(L,h,J,flags.ext);

%% ----- step 3 : Run computation 
c = comp_fwt(f,h.filts,J,h.a,Lc,flags.ext);

%% ----- FINALIZE : Reshape back 
permutedsizeAlt = size(c);
c=assert_sigreshape_post(c,dim,permutedsizeAlt,order);

%END FWT
 
 


