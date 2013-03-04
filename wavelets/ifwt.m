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
%   by the subsampling operation.
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
[flags,kv,Ls,dim]=ltfatarghelper({'Ls','dim'},definput,varargin);

% Initialize the wavelet filters structure
g = fwtinit(g,'syn');


Lc = [];
% If cell array was passed change the format.
if iscell(c)
   [c,Lc] =  wavcell2pack(c,dim);
end

%% ----- step 0 : Change c to correct shape according to the dim. 
[c,~,Lcsum,W,dim,permutedsize,order]=assert_sigreshape_pre(c,[],dim,upper(mfilename));


%% ----- step 1 : Determine input data length.
L = fwtlength(Ls,g,J,flags.ext);

%% ----- step 2 : Determine number of ceoefficients in each subband
if isempty(Lc)
   Lc = fwtclength(L,g,J,flags.ext);
end

%% ----- step 3 : Run computation 
f = comp_ifwt(c,g.filts,J,g.a,Lc,Ls,flags.ext);

%% ----- FINALIZE: Reshape back according to the dim.
permutedsizeAlt = size(f);
f=assert_sigreshape_post(f,dim,permutedsizeAlt,order);

%END IFWT


