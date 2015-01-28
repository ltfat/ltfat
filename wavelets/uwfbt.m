function [c,info]=uwfbt(f,wt,varargin)
%UWFBT   Undecimated Wavelet FilterBank Tree
%   Usage:  c=uwfbt(f,wt);
%           [c,info]=uwfbt(...);
%
%   Input parameters:
%         f   : Input data.
%         wt  : Wavelet Filterbank tree
%
%   Output parameters:
%         c     : Coefficients stored in $L \times M$ matrix.
%
%   `uwfbt(f,wt)` computes redundant time (or shift) invariant 
%   representation of the input signal *f* using the filterbank tree 
%   definition in *wt* and using the "a-trous" algorithm. 
%   Number of columns in *c* (*M*) is defined by the total number of 
%   outputs of nodes of the tree.
%
%   `[c,info]=uwfbt(f,wt)` additionally returns struct. `info` containing
%   the transform parameters. It can be conviniently used for the inverse
%   transform |iuwfbt| e.g. `fhat = iuwfbt(c,info)`. It is also required 
%   by the |plotwavelets| function.
%
%   If *f* is a matrix, the transformation is applied to each of *W* columns
%   and the coefficients in *c* are stacked along the third dimension.
%
%   Please see help for |wfbt| description of possible formats of *wt* and
%   description of frequency and natural ordering of the coefficient subbands.
%
%   Filter scaling
%   --------------
%
%   When compared to |wfbt|, the subbands produced by |uwfbt| are
%   gradually more and more redundant with increasing depth in the tree.
%   This results in energy grow of the coefficients. There are 3 flags
%   defining filter scaling:
%
%      'sqrt'
%               Each filter is scaled by `1/sqrt(a)`, there *a* is the hop
%               factor associated with it. If the original filterbank is
%               orthonormal, the overall undecimated transform is a tight
%               frame.
%               This is the default.
%
%      'noscale'
%               Uses filters without scaling.
%
%      'scale'
%               Each filter is scaled by `1/a`.
%
%   If 'noscale' is used, 'scale' has to be used in |iuwfbt| (and vice
%   versa) in order to obtain a perfect reconstruction.
%
%   Examples:
%   ---------
%
%   A simple example of calling the |uwfbt| function using the "full decomposition" wavelet tree:::
%
%     f = greasy;
%     J = 8;
%     [c,info] = uwfbt(f,{'sym10',J,'full'});
%     plotwavelets(c,info,16000,'dynrange',90);
%
%   See also: iuwfbt, wfbtinit

% AUTHOR: Zdenek Prusa

complainif_notenoughargs(nargin,2,'UWFBT');

definput.import = {'wfbtcommon','uwfbtcommon'};
flags=ltfatarghelper({},definput,varargin);

% Initialize the wavelet tree structure
wt = wfbtinit(wt,flags.forder);

%% ----- step 1 : Verify f and determine its length -------
[f,Ls]=comp_sigreshape_pre(f,upper(mfilename),0);
if(Ls<2)
   error('%s: Input signal seems not to be a vector of length > 1.',upper(mfilename));
end

%% ----- step 2 : Prepare input parameters
[nodesBF, rangeLoc, rangeOut] = treeBFranges(wt);
nodesUps = nodesFiltUps(nodesBF,wt);
%% ----- step 3 : Run computation
c = comp_uwfbt(f,wt.nodes(nodesBF),nodesUps,rangeLoc,rangeOut,flags.scaling);

%% ----- Optional : Fill the info struct. -----
if nargout>1
   info.fname = 'uwfbt';
   info.wt = wt;
   info.fOrder = flags.forder;
   info.isPacked = 0;
   info.scaling = flags.scaling;
end
