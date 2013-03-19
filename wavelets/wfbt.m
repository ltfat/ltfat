function [c,info]=wfbt(f,wt,varargin)
%WFBT   Wavelet FilterBank Tree
%   Usage:  c=wfbt(f,wt);
%           c=wfbt(f,wt,...);
%
%   Input parameters:
%         f   : Input data.
%         wt  : Wavelet Filterbank tree
%
%   Output parameters:
%         c   : Coefficients stored in a cell-array.
%
%   `c=wfbt(f,wt)` returns coefficients *c* obtained applying wavelet filterbank tree
%   defined by *wt* to the input data *f*. If *f* is a matrix, the transformation 
%   is applied to each of *W* columns. The *wt* parameter can be structure
%   obtained from the |wfbtinit| function and modified arbitrarily or it
%   can be cell-array, which is used as a parameter in the internal call of
%   the |wfbtinit| function.
%
%   The following flag groups are supported:
%
%         'per','zero','odd','even'
%                Type of the boundary handling.
%
%         'dwt','full'
%                Type of the tree to be used.
%
%         'freq','nat'
%                Frequency or natural order of the coefficient subbands.
%
%   Please see the help on |fwt| for a description of the boundary condition flags.
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
%


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
   info.Lc = cellfun(@(cEl) size(c,1),c);
   info.dim = dim;
   info.Ls = Ls;
   info.fOrder = flags.forder;
   info.isPacked = 0;
end


