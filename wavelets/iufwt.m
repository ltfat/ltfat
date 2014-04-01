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
%   levels of the basic synthesis filterbank defined by *g* using the "a-trous"
%   algorithm. Node that the same flag as in the `ufwt` function have to be used.
%
%   Please see the help on |ufwt| for a description of the parameters.
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


if nargin<2
   error('%s: Too few input parameters.',upper(mfilename));
end;

if(isstruct(par)&&isfield(par,'fname'))
   if nargin>2
      error('%s: Too many input parameters.',upper(mfilename));
   end
   w = fwtinit({'dual',par.wt});
   J = par.J;
else
   if nargin<3
      error('%s: Too few input parameters.',upper(mfilename));
   end;
   definput.keyvals.J = [];
   [~,~,J]=ltfatarghelper({'J'},definput,varargin);

   complain_notposint(J,'J');
   
   % Initialize the wavelet filters structure
   w = fwtinit(par);
end


%%  Run computation
f = comp_iufwt(c,w.g,J,w.a);







