function [c,info]=uwpfbt(f,wt,varargin)
%UWPFBT Undecimated Wavelet Packet FilterBank Tree
%   Usage:  c=uwpfbt(f,wt);
%           c=uwpfbt(f,wt,...);
%
%   Input parameters:
%         f   : Input data.
%         wt  : Wavelet Filterbank tree
%
%   Output parameters:
%         c   : Coefficients in a 3D matrix.
%
%   `c=uwpfbt(f,wt)` returns wavelet packet coefficients *c* obtained by
%   applying an undecimated wavelet filterbank tree defined by *wt* to the
%   input data *f*. If *f* is a matrix, the transformation is applied to 
%   each of column of the matrix.
%
%   Examples:
%   ---------
%   
%   A simple example of calling the |uwpfbt| function using the "full
%   decomposition" wavelet tree:::
% 
%     f = greasy;
%     J = 6;
%     [c,info] = uwpfbt(f,{'db10',J,'full'});
%     plotwavelets(c,info,44100,'dynrange',90);
%
%   See also: iuwpfbt, wfbtinit


if(nargin<2)
   error('%s: Too few input parameters.',upper(mfilename));  
end

definput.import = {'wfbtcommon'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Initialize the wavelet tree structure
wt = wfbtinit(wt,flags.forder,'ana');
    
%% ----- step 1 : Verify f and determine its length -------
[f,Ls]=comp_sigreshape_pre(f,upper(mfilename),0);
if(Ls<2)
   error('%s: Input signal seems not to be a vector of length > 1.',upper(mfilename));  
end

%% ----- step 3 : Run computation
wtPath = nodesBForder(wt);
nodesUps = nodeFiltUps(wtPath,wt);
rangeLoc = rangeInLocalOutputs(wtPath,wt);
c = comp_uwpfbt(f,wt.nodes(wtPath),rangeLoc,nodesUps);

%% ----- Optional : Fill the info struct. -----
if nargout>1
   info.fname = 'uwpfbt';
   info.wt = wt;
   info.fOrder = flags.forder;
   info.isPacked = 0;
end