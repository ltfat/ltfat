function f = iufwt(c,par,varargin)
%IUFWT   Inverse Undecimated Fast Wavelet Transform
%   Usage:  f = iufwt(c,info)
%           f = iufwt(c,w,J);
%
%   Input parameters:
%         c      : Coefficients stored in $L \times J+1$ matrix.
%         info,w : Transform parameters struct/Wavelet filters definition.
%         J      : Number of filterbank iterations.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   `f = iufwt(c,info)` reconstructs signal *f* from the wavelet
%   coefficients *c* using parameters from `info` struct. both returned by
%   |ufwt| function.
%
%   `f = iufwt(c,w,J)` reconstructs signal *f* from the wavelet
%   coefficients *c* using the wavelet filterbank consisting of the *J*
%   levels of the basic synthesis filterbank defined by *w* using the "a-trous"
%   algorithm. Node that the same flag as in the `ufwt` function have to be used.
%
%   Please see the help on |ufwt| for a description of the parameters.
%
%   Filter scaling
%   --------------
%
%   As in |ufwt|, 3 flags defining scaling of filters are recognized:
%
%      'sqrt'
%               Each filter is scaled by `1/sqrt(a)`, there *a* is the hop
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
%   If 'noscale' is used, 'scale' must have been used in |ufwt| (and vice
%   versa) in order to obtain a perfect reconstruction.
%
%   Examples:
%   ---------
%
%   A simple example showing perfect reconstruction:::
%
%     f = gspi;
%     J = 8;
%     c = ufwt(f,'db8',J);
%     fhat = iufwt(c,'db8',J);
%     % The following should give (almost) zero
%     norm(f-fhat)
%
%   See also:  ufwt, plotwavelets
%
%   References: ma98

% AUTHOR: Zdenek Prusa

complainif_notenoughargs(nargin,2,'IUFWT');

if(isstruct(par)&&isfield(par,'fname'))
   complainif_toomanyargs(nargin,2,'IUFWT');
   
   if ~strcmpi(par.fname,'ufwt')
      error(['%s: Wrong func name in info struct. ',...
             ' The info parameter was created by %s.'],...
             upper(mfilename),par.fname);
   end

   % Ensure we are using the correct wavelet filters
   w = fwtinit({'dual',par.wt});
   J = par.J;
   scaling = par.scaling;
   % Use the "oposite" scaling
   if strcmp(scaling,'scale')
       scaling = 'noscale';
   elseif strcmp(scaling,'noscale')
       scaling = 'scale';
   end
else
   complainif_notenoughargs(nargin,3,'IUFWT');

   definput.keyvals.J = [];
   definput.import = {'uwfbtcommon'};
   [flags, ~, J]=ltfatarghelper({'J'},definput,varargin);
   complainif_notposint(J,'J');

   % Initialize the wavelet filters
   % It is up to user to use the correct ones.
   w = fwtinit(par);
   scaling = flags.scaling;
end


%%  Run computation
f = comp_iufwt(c,w.g,w.a,J,scaling);







