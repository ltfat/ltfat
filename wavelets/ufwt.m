function [c,info] = ufwt(f,h,J,varargin)
%UFWT  Undecimated Fast Wavelet Transform 
%   Usage:  c = ufwt(f,h,J);
%
%   Input parameters:
%         f     : Input data.
%         h     : Analysis Wavelet Filterbank. 
%         J     : Number of filterbank iterations.
%
%   Output parameters:
%         c      : Coefficients.
%
%   `c=ufwt(f,h,J)` computes redundant time (or shift) invariant
%   wavelet representation *c* of the input signal *f* using the "a-trous"
%   algorithm.
%
%   The coefficents *c* are so called Undecimated Discrete Wavelet transform
%   of the input signal *f*, if *h* defines two-channel wavelet filterbank.
%   Other names for this version of the wavelet transform are: the
%   time-invariant wavelet transform, the stationary wavelet transform, 
%   maximal overlap discrete wavelet transform or even the "continuous" 
%   wavelet transform (as the time step is one sample).
%   However, the function accepts any number channels (further referred to as $N$)
%   in the basic wavelet filterbank. 
%   The redundancy of the representanion is $J*(N-1)+1$.
%
%   For all accepted formats of the parameter *h* see the |fwt| function.
%
%   For one-dimensional input *f*, the coefficients *c* are stored as collumns of
%   a matrix. The collumns are ordered with inceasing central frequency of the
%   corresponding effective filter frequency response or equivalently with decreasing wavelet scale.
%
%   If the input *f* is matrix with *W* columns, the transform is applied to each collumn
%   and the outputs are stacked along third dimension.
%
%
%   Boundary handling:
%   ------------------
%
%   `c=ufwt(f,h,J)` uses periodic boundary extension. The extensions are 
%   done at each level of the transform internally rather than doing the 
%   prior explicit padding.
%   
%
%   Examples:
%   ---------
%   
%   A simple example of calling the |ufwt| function:::
% 
%     f = gspi;
%     J = 8;
%     [c,info] = ufwt(f,'db8',J);
%     %plotwavelets(c,info,44100,'dynrange',90);
%
%   See also: iufwt, fwtinit
%
%   References: holschneider1989real

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
h = fwtinit(h,'ana');

%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
[f,Ls]=comp_sigreshape_pre(f,upper(mfilename),0);
if(Ls<2)
   error('%s: Input signal seems not to be a vector of length > 1.',upper(mfilename));  
end
 
%% ----- step 2 : Run computation
 c = comp_ufwt(f,h.filts,J,h.a);
 
%% ----- Optionally : Fill info struct ----
if nargout>1
   info.fname = 'ufwt';
   info.fwtstruct = h;
   info.J = J;
end


