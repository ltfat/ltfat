function [c,info]=uwfbt(f,wt,varargin)
%UWFBT   Undecimated Wavelet FilterBank Tree
%   Usage:  c=uwfbt(f,wt);
%           c=uwfbt(f,wt,...);
%
%   Input parameters:
%         f   : Input data.
%         wt  : Wavelet Filterbank tree
%
%   Output parameters:
%         c   : Coefficients stored in a cell-array.
%
%   `c=uwfbt(f,wt)` returns coefficients *c* obtained by applying wavelet filterbank tree
%   defined by *wt* to the input data *f*. If *f* is a matrix, the transformation 
%   is applied to each of *W* columns. 
%
%   This routine accepts the same flags as |fwt| for specifying the boundary
%   conditions. Please see the help on |fwt| for a description of the flags.
%
%   Examples:
%   ---------
%   
%   A simple example of calling the |wfbt| function using the "full decomposition" wavelet tree:::
% 
%     f = greasy;
%     J = 8;
%     [c,info] = uwfbt(f,{'sym10',J,'full'});
%     plotwavelets(c,info,16000,'dynrange',90);
%
%   See also: iuwfbt, wfbtinit

if(nargin<2)
   error('%s: Too few input parameters.',upper(mfilename));  
end

definput.import = {'wfbtcommon'};
flags=ltfatarghelper({},definput,varargin);

% Initialize the wavelet tree structure
wt = wfbtinit(wt,flags.forder,'ana');

%% ----- step 1 : Verify f and determine its length -------
[f,Ls]=comp_sigreshape_pre(f,upper(mfilename),0);
if(Ls<2)
   error('%s: Input signal seems not to be a vector of length > 1.',upper(mfilename));  
end

%% ----- step 2 : Prepare input parameters
wtPath = nodesBForder(wt);
nodesUps = nodeFiltUps(wtPath,wt);
rangeLoc = rangeInLocalOutputs(wtPath,wt);
rangeOut = rangeInOutputs(wtPath,wt);
%% ----- step 3 : Run computation
c = comp_uwfbt(f,wt.nodes(wtPath),nodesUps,rangeLoc,rangeOut);

%% ----- Optional : Fill the info struct. -----
if nargout>1
   info.fname = 'uwfbt';
   info.wt = wt;
   info.fOrder = flags.forder;
   info.isPacked = 0;
end
