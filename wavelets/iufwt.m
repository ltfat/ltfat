function f = iufwt(c,par,varargin)
%IUFWT   Inverse Undecimated Fast Wavelet Transform 
%   Usage:  f = iufwt(c,info)
%           f = iufwt(c,g,J);   
%
%   Input parameters:
%         c     : Coefficients stored in a cell-array.
%         g     : Synthesis wavelet filters.
%         J     : Number of filterbank iterations.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   `f = iufwt(c,g,J)` reconstruct signal *f* from the wavelet
%   coefficients *c* using the wavelet filterbank consisting of the *J* levels
%   of the basic synthesis filterbank defined by *g* using the "a-trous"
%   algorithm.
%
%   Node that the same flag as in the `ufwt` function have to be used.
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
%   See also:  ufwt, fwtinit
%
%   References: ma98


if nargin<2
   error('%s: Too few input parameters.',upper(mfilename));
end;

if(isstruct(par)&&isfield(par,'fname'))
   if nargin>2
      error('%s: Too many input parameters.',upper(mfilename));
   end
   g = fwtinit(par.fwtstruct,'syn');
   J = par.J;
else
   if nargin<3
      error('%s: Too few input parameters.',upper(mfilename));
   end;
   definput.keyvals.J = [];
   [~,~,J]=ltfatarghelper({'J'},definput,varargin);

   if ~isnumeric(J) || ~isscalar(J)
     error('%s: "J" must be a scalar.',upper(mfilename));
   end;

   if(J<1 && rem(a,1)~=0)
      error('%s: J must be a positive integer.',upper(mfilename)); 
   end
   
   % Initialize the wavelet filters structure
   g = fwtinit(par,'syn');
end


%%  Run computation
f = comp_iufwt(c,g.g,J,g.a);







