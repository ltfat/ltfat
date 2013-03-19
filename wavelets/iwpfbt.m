function f=iwpfbt(c,par,varargin)
%IWPFBT   Inverse Wavelet Packet Filterbank Tree
%   Usage:  f=iwpfbt(c,info);
%           f=iwpfbt(c,wt,Ls,...);
%
%   Input parameters:
%         c     : Coefficients stored in a cell-array.
%         wt    : Wavelet Filterbank tree
%         Ls    : Length of the reconstructed signal.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   `f=iwpfbt(c,wt)` reconstructs signal *f* from coefficients *c* using the
%   wavelet filterbank tree *wt*. 
%
%   The following flag groups are supported:
%
%         'per','zpd','sym','symw','asym','asymw','ppd','sp0'
%                Type of the boundary handling.
%
%         'full','dwt'
%                Type of the tree to be used.
%
%         'freq','nat'
%                Frequency or natural order of the coefficient subbands.
%
%   Please see the help on |fwt| for a description of the flags.
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
%   See also: wfbt, wfbtinit
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
   wt = wfbtinit(par.wt,par.fOrder,'syn');
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
   wt = wfbtinit(par,flags.forder,'syn');
end

wtPath = nodesBForder(wt,'rev');
[pOutIdxs,chOutIdxs] = rangeWpBF(wt,'rev');
f = comp_iwpfbt(c,wt.nodes(wtPath),pOutIdxs,chOutIdxs,Ls,ext);
