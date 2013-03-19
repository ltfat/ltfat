function f = ifwt(c,par,varargin)
%IFWT   Inverse Fast Wavelet Transform 
%   Usage:  f = ifwt(c,info)
%           f = ifwt(c,g,J,Ls,...)
%
%   Input parameters:
%         c     : Wavelet coefficients.
%         g     : Synthesis wavelet filters.
%         J     : Number of filterbank iterations.
%         Ls    : Length of the reconstructed signal.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   `f = ifwt(c,g,J)` reconstructs signal *f* from the wavelet coefficients
%   *c* using *J*-iteration synthesis filter bank build from the basic synthesis
%   filterbank defined by *g*. The fast wavelet transform algorithm 
%   (Mallat's algorithm) is employed. The format of *c* can be either
%   packed, as returned by the |fwt| function or cell-array as returned by
%   |wavpack2cell| function.
%
%   The *Ls* parameter is mandatory due to the ambiguity of lengths introduced
%   by the subsampling operation and by boundary treatment methods.
%
%   `ifwt` supports the same boundary conditions as |fwt|. Note that the
%   same flag as in the |fwt| function have to be used, otherwise
%   perfect reconstruction cannot be obtained.
%
%   Please see the help on |fwt| for a description of the parameters.
%
%   Examples:
%   ---------
%   
%   A simple example showing perfect reconstruction:::
% 
%     f = gspi;
%     J = 8;
%     c = fwt(f,'db8',J);
%     fhat = ifwt(c,'db8',J,length(f));
%     % The following should give (almost) zero
%     norm(f-fhat)
%   
%   See also:  fwt, fwtinit
%
%   References: ma98

if nargin<2
   error('%s: Too few input parameters.',upper(mfilename));
end;

if  ~(iscell(c) || isnumeric(c))
  error('%s: Unrecognized coefficient format.',upper(mfilename));
end

if(isstruct(par)&&isfield(par,'fname'))
   if nargin>2
      error('%s: Too many input parameters.',upper(mfilename));
   end
   % process info struct
   g = fwtinit(par.fwtstruct,'syn');
   J = par.J;
   Lc = par.Lc;
   Ls = par.Ls;
   dim = par.dim;
   ext = par.ext;
else
   if nargin<4
      error('%s: Too few input parameters.',upper(mfilename));
   end;
   
   %% PARSE INPUT
   definput.import = {'fwt'};
   definput.keyvals.dim = [];
   definput.keyvals.Ls = [];
   definput.keyvals.J = [];
   [flags,~,J,Ls,dim]=ltfatarghelper({'J','Ls','dim'},definput,varargin);
 
   ext = flags.ext;
   %If dim is not specified use first non-singleton dimension.
   if(isempty(dim))
       dim=find(size(c)>1,1);
   end
   
   % Initialize the wavelet filters structure
   g = fwtinit(par,'syn');

   %% ----- Determine input data length.
   L = fwtlength(Ls,g,J,ext);

   %% ----- Determine number of ceoefficients in each subband
   Lc = fwtclength(L,g,J,ext);
end

   %% ----- Change c to correct shape according to the dim.
   if(isnumeric(c))
      %Change format
      c =  wavpack2cell(c,Lc,dim);
   elseif(iscell(c))
      %Just check *Lc*
      if(~isequal(Lc,cellfun(@(x) size(x,1),c)))
         error('%s: Coefficient subband lengths does not comply with parameter *Ls*.',upper(mfilename));
      end
   else
      error('%s: Unrecognized coefficient format.',upper(mfilename)); 
   end;

   %% ----- Run computation 
   f = comp_ifwt(c,g.filts,J,g.a,Ls,ext);

   %% ----- FINALIZE: Reshape back according to the dim.
   if(dim==2)
       f = f.';
   end
%END IFWT
