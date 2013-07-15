function f=iwpfbt(c,par,varargin)
%IWPFBT   Inverse Wavelet Packet Filterbank Tree
%   Usage:  f=iwpfbt(c,info);
%           f=iwpfbt(c,wt,Ls);
%
%   Input parameters:
%         c       : Coefficients stored in a cell-array.
%         info,wt : Transform parameters struct/Wavelet Filterbank tree.
%         Ls      : Length of the reconstructed signal.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   `f = iwpfbt(c,info)` reconstructs signal *f* from the coefficients *c* 
%   using parameters from `info` struct. both returned by |wfbt| function.
%
%   `f = iwpfbt(c,wt,Ls)` reconstructs signal *f* from the coefficients *c*
%   using filter bank tree defined by *wt*. Plese see |wfbt| function for
%   possible formats of *wt*. The *Ls* parameter is mandatory due to the 
%   ambiguity of reconstruction lengths introduced by the subsampling 
%   operation and by boundary treatment methods. Note that the same flag as
%   in the |wfbt| function have to be used, otherwise perfect reconstruction
%   cannot be obtained. Please see help for |wfbt| for description of the
%   flags.
%
%   Examples:
%   ---------
%   
%   A simple example showing perfect reconstruction using the "full decomposition" wavelet tree:::
% 
%     f = gspi;
%     J = 7;
%     wtdef = {'db10',J,'full'};
%     c = wpfbt(f,wtdef);
%     fhat = iwpfbt(c,wtdef,length(f));
%     % The following should give (almost) zero
%     norm(f-fhat)
%
%   See also: wpfbt, wfbtinit
%
if nargin<2
   error('%s: Too few input parameters.',upper(mfilename));
end;

if(~iscell(c))
    error('%s: Unrecognized coefficient format.',upper(mfilename));
end

if(isstruct(par)&&isfield(par,'fname'))
   if nargin>2
      error('%s: Too many input parameters.',upper(mfilename));
   end
   wt = wfbtinit({'dual',par.wt},par.fOrder);
   Ls = par.Ls;
   ext = par.ext;
else
   if nargin<3
      error('%s: Too few input parameters.',upper(mfilename));
   end
   %% PARSE INPUT
   definput.keyvals.Ls=[];    
   definput.import = {'fwt','wfbtcommon'};
   [flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);
   ext = flags.ext;
   % Initialize the wavelet tree structure
   wt = wfbtinit(par,flags.forder);
end

wtPath = nodesBForder(wt,'rev');
[pOutIdxs,chOutIdxs] = rangeWpBF(wt,'rev');
f = comp_iwpfbt(c,wt.nodes(wtPath),pOutIdxs,chOutIdxs,Ls,ext);
