function f=iuwpfbt(c,wt,varargin)
%IWPFBT Inverse Undecimated Wavelet Packet Filterbank Tree
%   Usage:  f=iwpfbt(c,wt);
%           f=iwpfbt(c,wt,Ls,...);
%
%   Input parameters:
%         c     : Coefficients stored in a cell-array.
%         wt    : Wavelet Filterbank tree
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   `f=iuwpfbt(c,wt)` reconstructs signal *f* from coefficients *c* using the
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
%     fhat = iuwpfbt(c,wt,length(f));
%     % The following should give (almost) zero
%     norm(f-fhat)
%
%   See also: wfbt, wfbtinit
%
if nargin<2
   error('%s: Too few input parameters.',upper(mfilename));
end;


%% PARSE INPUT
definput.import = {'wfbtcommon'};
if(~isnumeric(c))
    error('%s: Unrecognized coefficient format.',upper(mfilename));
end

[flags,kv]=ltfatarghelper({},definput,varargin);

% Initialize the wavelet tree structure
wt = wfbtinit(wt,flags.forder,'syn');


wtPath = nodesBForder(wt,'rev');
[pOutIdxs,chOutIdxs] = rangeWpBF(wt,'rev');
nodesUps = nodeFiltUps(wtPath,wt);
f = comp_iuwpfbt(c,wt.nodes(wtPath),nodesUps,pOutIdxs,chOutIdxs);
