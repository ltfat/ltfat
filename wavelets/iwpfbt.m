function f=iwpfbt(c,wt,Ls,varargin)
%IWPFBT   Inverse Wavelet Packet Filterbank Tree
%   Usage:  f=iwpfbt(c,wt,Ls);
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
%   Please see the help on |fwt|_ for a description of the flags.
%
%   Examples:
%   ---------
%   
%   A simple example showing perfect reconstruction using the "full decomposition" wavelet tree:::
% 
%     f = gspi;
%     J = 7;
%     wt = wfbtinit({{'db',10},J},'full');
%     c = wpfbt(f,wt);
%     fhat = iwpfbt(c,wt,length(f));
%     % The following should give (almost) zero
%     norm(f-fhat)
%
%   See also: wfbt, wfbtinit
%
if nargin<3
   error('%s: Too few input parameters.',upper(mfilename));
end;


%% PARSE INPUT
definput.keyvals.Ls=[];    
definput.import = {'fwt','wfbtcommon'};

if(~iscell(c))
    error('%s: Unrecognized coefficient format.',upper(mfilename));
end

[flags,kv]=ltfatarghelper({},definput,varargin);

% Initialize the wavelet tree structure
wt = wfbtinit(wt,flags.treetype,'syn');

f = comp_iwpfbt(c,wt,Ls,'dec',flags.ext);
