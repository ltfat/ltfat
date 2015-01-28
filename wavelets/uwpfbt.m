function [c,info]=uwpfbt(f,wt,varargin)
%UWPFBT Undecimated Wavelet Packet FilterBank Tree
%   Usage:  c=uwpfbt(f,wt);
%           [c,info]=uwpfbt(...);
%
%   Input parameters:
%         f   : Input data.
%         wt  : Wavelet Filterbank tree
%
%   Output parameters:
%         c   : Coefficients in a $L \times M$ matrix.
%
%   `c=uwpfbt(f,wt)` returns coefficients *c* obtained by applying the
%   undecimated wavelet filterbank tree defined by *wt* to the input data
%   *f* using the "a-trous" algorithm. Number of columns in *c* (*M*) is
%   defined by the total number of outputs of each node. The outputs `c(:,jj)`
%   are ordered in the breadth-first node order manner.
%
%   `[c,info]=uwpfbt(f,wt)` additionally returns struct. `info` containing 
%   the transform parameters. It can be conviniently used for the inverse 
%   transform |iuwpfbt| e.g. `fhat = iuwpfbt(c,info)`. It is also required
%   by the |plotwavelets| function.
%
%   If *f* is a matrix, the transformation is applied to each of *W* columns
%   and the coefficients in *c* are stacked along the third dimension.
%
%   Please see help for |wfbt| description of possible formats of *wt*.
%
%   Scaling of intermediate outputs:
%   --------------------------------
%
%   The following flags control scaling of intermediate outputs and
%   therefore the energy relations between coefficient subbands. An 
%   intermediate output is an output of a node which is further used as an
%   input to a descendant node.
%
%      'intsqrt'
%               Each intermediate output is scaled by `1/sqrt(2)`.
%               If the filterbank in each node is orthonormal, the overall
%               undecimated transform is a tight frame.
%               This is the default.
%
%      'intnoscale'
%               No scaling of intermediate results is used. This is
%               necessaty for the |wpbest| function to correctly work with
%               the cost measures.
%
%      'intscale'
%               Each intermediate output is scaled by `1/2`.
%
%   If 'intnoscale' is used, 'intscale' must be used in |iuwpfbt| (and vice
%   versa) in order to obtain a perfect reconstruction.
%
%   Scaling of filters:
%   -------------------
%
%   When compared to |wpfbt|, the subbands produced by |uwpfbt| are
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
%   If 'noscale' is used, 'scale' must be used in |iuwpfbt| (and vice
%   versa) in order to obtain a perfect reconstruction.
%
%   Examples:
%   ---------
%
%   A simple example of calling the |uwpfbt| function using the "full
%   decomposition" wavelet tree:::
%
%     [f,fs] = greasy;
%     J = 6;
%     [c,info] = uwpfbt(f,{'db10',J,'full'});
%     plotwavelets(c,info,fs,'dynrange',90);
%
%   See also: iuwpfbt, wfbtinit

complainif_notenoughargs(nargin,2,'UWPFBT');

definput.import = {'wfbtcommon','uwfbtcommon'};
definput.flags.interscaling = {'intsqrt', 'intscale', 'intnoscale'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Initialize the wavelet tree structure
wt = wfbtinit(wt,flags.forder);

%% ----- step 1 : Verify f and determine its length -------
[f,Ls]=comp_sigreshape_pre(f,upper(mfilename),0);
if(Ls<2)
   error('%s: Input signal seems not to be a vector of length > 1.',...
         upper(mfilename));
end

%% ----- step 3 : Run computation
wtPath = nodeBForder(0,wt);
nodesUps = nodesFiltUps(wtPath,wt);
rangeLoc = nodesLocOutRange(wtPath,wt);
c = comp_uwpfbt(f,wt.nodes(wtPath),rangeLoc,nodesUps,flags.scaling,...
                flags.interscaling);

%% ----- Optional : Fill the info struct. -----
if nargout>1
   info.fname = 'uwpfbt';
   info.wt = wt;
   info.fOrder = flags.forder;
   info.isPacked = 0;
   info.interscaling = flags.interscaling;
   info.scaling = flags.scaling;
end
