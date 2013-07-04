function [c,info] = fwt(f,h,J,varargin)
%FWT   Fast Wavelet Transform 
%   Usage:  c = fwt(f,h,J);
%           c = fwt(f,h,J,'dim',dim,'ext',ext);
%           [c,info] = fwt(...);
%
%   Input parameters:
%         f     : Input data.
%         h     : Analysis Wavelet Filterbank. 
%         J     : Number of filterbank iterations.
%
%   Output parameters:
%         c      : Coefficient vector.
%         info   : Transform parameters struct.
%
%   `c=fwt(f,h,J)` returns wavelet coefficients *c* of the input signal *f*
%   using *J* iterations of the basic wavelet filterbank defined by *h*
%   i.e. the fast wavelet transform algorithm (Mallat's algorithm) is used.
%   In addition, the function returns struct. `info` containing transform
%   parameters. It can be conviniently used for the inverse transform 
%   |ifwt| e.g. `fhat = ifwt(c,info)`. It is also required by the 
%   |plotwavelets| function.
%   
%   The coefficents *c* are Discrete Wavelet transform of the input signal
%   *f*, if *h* defines two-channel wavelet filterbank. The function can 
%   apply the Mallat's algorithm using basic filter banks with any number 
%   of the channels. In such case, the transform have a different name.    
%
%   The basic analysis wavelet filterbank $h$ can be passed in several 
%   formats. The formats are the same as for the |fwtinit| function.
%
%   1) One option is passing a cell array whose first element is the
%      name of the function defining the basic wavelet filters (`wfilt_`
%      prefix) and the other elements are the parameters of the
%      function. e.g. `{'db',10}` calls `wfilt_db(10)` internally.
%
%   2) Character string as concatenation of the name of the wavelet
%      filters defining function (as above) and the numeric parameters
%      delimited by ':' character, e.g. 'db10' has the same effect as above,
%      'spline4:4' calls `wfilt_spline(4,4)` internally.
%
%   3) The third possible format of $h$ is to pass cell array of one
%      dimensional numerical vectors directly defining the wavelet filter
%      impulse responses.  In this case, outputs of the filters are
%      subsampled by a factor equal to the number of the filters. For
%      creating completely custom filterbanks use the |fwtinit| function.
%
%   4) The fourth option is to pass a structure obtained from the
%      |fwtinit| function.
%   
%   If *f* is row/collumn vector of length *L*, the subbands *c* are stored
%   in a single collumn/row in a consecutive order with respect to the 
%   inceasing central frequency of the corresponding effective filter 
%   frequency responses or equivalently with decreasing wavelet scale. The 
%   lengths of subbands are stored in *info.Lc* so the subbands can be easily
%   extracted using |wavpack2cell|.
%
%   If the input *f* is a matrix, the transform is applied to each column 
%   if `dim==1` (default) and `Ls=size(f,1)` and `W=size(f,2)`. If `dim==2`
%   the transform is applied to each row `W=size(f,1)` and `Ls=size(f,2)`.
%   The output is then a matrix and the input orientation is preserved in 
%   the orientation of the output coefficients.
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
%   `fwt(f,h,J,'ext',ext)` with `ext` other than `'per'` computes a slightly
%   redundant wavelet representation of the input signal *f* with the chosen
%   boundary extension *ext*. The redundancy (expansivity) of the
%   represenation is the price to pay for using general filterbank and
%   custom boundary treatment.  The extensions are done at each level of the
%   transform internally rather than doing the prior explicit padding.
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
%   function |ifwt| if the `info` struct is not used.
%
%   Examples:
%   ---------
%   
%   A simple example of calling the |fwt| function:::
% 
%     f = greasy;
%     J = 10;
%     [c,info] = fwt(f,'db8',J);
%     plotwavelets(c,info,44100,'dynrange',90);
%
%   See also: ifwt, fwtinit, wavpack2cell, plotwavelets
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
definput.flags.cfmt = {'pack','cell'};
[flags,~,dim]=ltfatarghelper({'dim'},definput,varargin);


%% ----- step 1 : Verify f and determine its length -------
[f,~,Ls,~,dim]=assert_sigreshape_pre(f,[],dim,upper(mfilename));

%% ----- step 2 : Determine number of ceoefficients in each subband *Lc*
%  and next legal input data length *L*.
[Lc, L] = fwtclength(Ls,h,J,flags.ext);

% Pad with zeros if the safe length L differ from the Ls.
if(Ls~=L)
   f=postpad(f,L); 
end


%% ----- step 3 : Run computation. 
c = comp_fwt(f,h.filts,J,h.a,Lc,flags.ext);

%% ----- FINALIZE: Change format of coefficients.
if flags.do_pack
   c = wavcell2pack(c,dim);
end

%% ----- FILL INFO STRUCT ----------------------
if nargout>1
   info.fname = 'fwt';
   info.fwtstruct = h;
   info.ext = flags.ext;
   info.Lc = Lc;
   info.J = J;
   info.dim = dim;
   info.Ls = Ls;
   info.isPacked = flags.do_pack;
end
%END FWT
