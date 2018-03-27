function [c,info] = ufwt(f,w,J,varargin)
%UFWT  Undecimated Fast Wavelet Transform
%   Usage:  c = ufwt(f,w,J);
%           [c,info] = ufwt(...);
%
%   Input parameters:
%         f     : Input data.
%         w     : Wavelet Filterbank.
%         J     : Number of filterbank iterations.
%
%   Output parameters:
%         c     : Coefficients stored in $L \times (J+1)$ matrix.
%         info  : Transform paramaters struct.
%
%   `ufwt(f,w,J)` computes redundant time (or shift) invariant
%   wavelet representation of the input signal *f* using wavelet filters
%   defined by *w* in the "a-trous"  algorithm. 
%
%   For all accepted formats of the parameter *w* see the |fwtinit| function.
%
%   `[c,info]=ufwt(f,w,J)` additionally returns the `info` struct. 
%   containing the transform parameters. It can be conviniently used for 
%   the inverse transform |iufwt| e.g. `fhat = iufwt(c,info)`. It is also 
%   required by the |plotwavelets| function.
%
%   The coefficents *c* are so called undecimated Discrete Wavelet transform
%   of the input signal *f*, if *w* defines two-channel wavelet filterbank.
%   Other names for this version of the wavelet transform are: the
%   time-invariant wavelet transform, the stationary wavelet transform,
%   maximal overlap discrete wavelet transform or even the "continuous"
%   wavelet transform (as the time step is one sample). However, the
%   function accepts any number filters (referred to as $M$) in the basic
%   wavelet filterbank and the number of columns of *c* is then $J(M-1)+1$.
%
%   For one-dimensional input *f* of length *L*, the coefficients *c* are
%   stored as columns of a matrix. The columns are ordered with inceasing
%   central frequency of the respective subbands.
%
%   If the input *f* is $L \times W$ matrix, the transform is applied
%   to each column and the outputs are stacked along third dimension in the
%   $L \times J(M-1)+1 \times W$ data cube.
%
%   Filter scaling
%   --------------
%
%   When compared to |fwt|, |ufwt| subbands are gradually more and more 
%   redundant with increasing level of the subband. If no scaling of the 
%   filters is introduced, the energy of subbands tends to grow with increasing
%   level.
%   There are 3 flags defining filter scaling:
%
%      'sqrt'
%               Each filter is scaled by `1/sqrt(a)`, where *a* is the hop
%               factor associated with it. If the original filterbank is
%               orthonormal, the overall undecimated transform is a tight
%               frame.
%               This is the default.
%
%      'noscale'
%               Uses filters without scaling.
%
%      'scale'
%               Each filter is scaled by `1/a`.
%
%   If 'noscale' is used, 'scale' has to be used in |iufwt| (and vice
%   versa) in order to obtain a perfect reconstruction.
%
%   Boundary handling:
%   ------------------
%
%   `c=ufwt(f,w,J)` uses periodic boundary extension. The extensions are
%   done internally at each level of the transform, rather than doing the
%   prior explicit padding.
%
%   Examples:
%   ---------
%
%   A simple example of calling the |ufwt| function:::
%
%     [f,fs] = greasy;
%     J = 8;
%     [c,info] = ufwt(f,'db8',J);
%     plotwavelets(c,info,fs,'dynrange',90);
%
%   See also: iufwt, plotwavelets
%
%   References: holschneider1989real

% AUTHOR: Zdenek Prusa

complainif_notenoughargs(nargin,3,'UFWT');
complainif_notposint(J,'J');

definput.import = {'uwfbtcommon'};
[flags]=ltfatarghelper({},definput,varargin);

% Initialize the wavelet filters structure
w = fwtinit(w);

%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
[f,Ls]=comp_sigreshape_pre(f,upper(mfilename),0);
if(Ls<2)
   error('%s: Input signal seems not to be a vector of length > 1.',upper(mfilename));
end

%% ----- step 2 : Run computation
c = comp_ufwt(f,w.h,w.a,J,flags.scaling);

%% ----- Optionally : Fill info struct ----
if nargout>1
   info.fname = 'ufwt';
   info.wt = w;
   info.J = J;
   info.scaling = flags.scaling;
end


