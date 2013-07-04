function [c,info]=wfbt(f,wt,varargin)
%WFBT   Wavelet FilterBank Tree
%   Usage:  c=wfbt(f,wt,...);
%           [c,info]=wfbt(...);
%
%   Input parameters:
%         f   : Input data.
%         wt  : Wavelet Filterbank tree definition.
%
%   Output parameters:
%         c    : Coefficients stored in a cell-array.
%         info : Transform parameters struct.
%
%   `c=wfbt(f,wt)` returns coefficients *c* obtained applying general wavelet 
%   filterbank tree defined by *wt* to the input data *f*. In addition, the
%   function returns struct. `info` containing transform parameters. It can
%   be conviniently used for the inverse transform |iwfbt| e.g. 
%   `fhat = iwfbt(c,info)`. It is also required by the |plotwavelets| function.
%   If *f* is a matrix, the transformation is applied to each of *W* columns. 
%
%   The *wt* parameter can have two formats:
%
%   1) Cell array containing 3 elements `{w,J,treetype}`, where `w` is
%      the basic wavelet filterbank definition as in |fwt| function, *J*
%      stands for the depth of the tree and the flag `treetype` defines 
%      the type of the tree to be used. Supported options are:
%
%      `'dwt'`  
%        DWT tree. Just the low-pass output is decomposed further.
%
%      `'full'`
%        Full decomposition tree. Each output is decomposed up to level *J*.
%
%   2) Structure returned by the |wfbtinit| function and possibly
%      modified by |wfbtput| and |wfbtremove|.
%
%   In addition, the following flag groups are supported:
%
%   `'per'`,`'zero'`,`'odd'`,`'even'`
%     Type of the boundary handling.
%
%   `'freq'`,`'nat'`
%     Frequency or natural order of the coefficient subbands.
%
%   Please see the help on |fwt| for a description of the boundary condition flags.
%
%   The output coefficients are stored in a cell array. Each element is a
%   2D matrix with *W* columns each of which respresents one output subband
%   of the w'th input channel.
%
%   Examples:
%   ---------
%   
%   A simple example of calling the |wfbt| function using the "full decomposition" wavelet tree:::
% 
%     f = gspi;
%     J = 7;
%     [c,info] = wfbt(f,{'sym10',J,'full'});
%     plotwavelets(c,info,44100,'dynrange',90);
%
%   See also: iwfbt, wfbtinit


if(nargin<2)
   error('%s: Too few input parameters.',upper(mfilename));  
end

definput.import = {'fwt','wfbtcommon'};
definput.keyvals.dim = [];
[flags,kv,dim]=ltfatarghelper({'dim'},definput,varargin);

% Initialize the wavelet tree structure
wt = wfbtinit(wt,flags.forder,'ana');
    
%% ----- step 1 : Verify f and determine its length -------
[f,~,Ls]=assert_sigreshape_pre(f,[],dim,upper(mfilename));

% Determine next legal input data length.
L = wfbtlength(Ls,wt,flags.ext);

% Pad with zeros if the safe length L differ from the Ls.
if(Ls~=L)
   f=postpad(f,L); 
end

%% ----- step 3 : Run computation
wtPath = nodesBForder(wt);
rangeLoc = rangeInLocalOutputs(wtPath,wt);
rangeOut = rangeInOutputs(wtPath,wt); % very slow
c = comp_wfbt(f,wt.nodes(wtPath),rangeLoc,rangeOut,flags.ext);

%% ----- Optionally : Fill info struct ----
if nargout>1
   info.fname = 'wfbt';
   info.wt = wt;
   info.ext = flags.ext;
   info.Lc = cellfun(@(cEl) size(cEl,1),c);
   info.dim = dim;
   info.Ls = Ls;
   info.fOrder = flags.forder;
   info.isPacked = 0;
end


