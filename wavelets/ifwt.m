function f = ifwt(c,g,J,varargin)
%IFWT   Inverse Fast Wavelet Transform 
%   Usage:  f = ifwt(c,g,J,Ls)
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
%   packed, as returned by the |fwt|_ function or cell-array as returned by
%   |wavpack2cell|_ function.
%
%   The *Ls* parameter is mandatory due to the ambiguity of lengths introduced
%   by the subsampling operation and by boundary treatment methods.
%
%   `ifwt` supports the same boundary conditions as |fwt|_. Note that the
%   same flag as in the |fwt|_ function have to be used, otherwise
%   perfect reconstruction cannot be obtained.
%
%   Please see the help on |fwt|_ for a description of the parameters.
%
%   Examples:
%   ---------
%   
%   A simple example showing perfect reconstruction:::
% 
%     f = gspi;
%     J = 8;
%     c = fwt(f,{'db',8},J);
%     fhat = ifwt(c,{'db',8},J,length(f));
%     % The following should give (almost) zero
%     norm(f-fhat)
%   
%   See also:  fwt, fwtinit
%
%   References: ma98

if nargin<4
   error('%s: Too few input parameters.',upper(mfilename));
end;


if  ~(iscell(c) || isnumeric(c))
  error('%s: Unrecognized coefficient format.',upper(mfilename));
end


%% PARSE INPUT
definput.import = {'fwt'};
definput.keyvals.dim = [];
definput.keyvals.Ls = [];
[flags,~,Ls,dim]=ltfatarghelper({'Ls','dim'},definput,varargin);

% Initialize the wavelet filters structure
g = fwtinit(g,'syn');

%If dim is not specified use first non-singleton dimension.
if(isempty(dim))
    dim=find(size(c)>1,1);
end

%% ----- step 1 : Determine input data length.
L = fwtlength(Ls,g,J,flags.ext);

%% ----- step 2 : Determine number of ceoefficients in each subband
Lc = fwtclength(L,g,J,flags.ext);

%% ----- step 4 : Change c to correct shape according to the dim.
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

%% ----- step 4 : Run computation 
f = comp_ifwt(c,g.filts,J,g.a,Ls,flags.ext);

%% ----- FINALIZE: Reshape back according to the dim.
if(dim==2)
    f = f.';
end

%END IFWT


