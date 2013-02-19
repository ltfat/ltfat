function c=wpfbt(f,wt,varargin)
%WPFBT   Wavelet Packet FilterBank Tree
%   Usage:  c=wpfbt(f,wt);
%           c=wpfbt(f,wt,...);
%
%   Input parameters:
%         f   : Input data.
%         wt  : Wavelet Filterbank tree
%
%   Output parameters:
%         c   : Coefficients stored in a cell-array.
%
%   `c=wpfbt(f,wt)` returns wavelet packet coefficients *c* obtained by
%   applying a wavelet filterbank tree defined by *wt* to the input data
%   *f*. If *f* is a matrix, the transformation is applied to each of column
%   of the matrix.
%
%   This routine supports the same boundary conditions flags as
%   |fwt|_. Please see the help on |fwt|_ for a description of the flags.
%
%   Examples:
%   ---------
%   
%   A simple example of calling the |wpfbt|_ function using the "full
%   decomposition" wavelet tree:::
% 
%     f = gspi;
%     J = 6;
%     c = wpfbt(f,{{'db',10},J},'full');
%     plotfwt(c,44100,90);
%
%   See also: iwfbt, wfbtinit


if(nargin<2)
   error('%s: Too few input parameters.',upper(mfilename));  
end

definput.import = {'fwt','wfbtcommon'};
%definput.flags.pack={'pack_null','pack'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Initialize the wavelet tree structure
wt = wfbtinit(wt,flags.treetype,flags.forder,'ana');

    
%% ----- step 1 : Verify f and determine its length -------
[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,upper(mfilename),0);
if(Ls<2)
   error('%s: Input signal seems not to be a vector of length > 1.',upper(mfilename));  
end

%% ----- step 2 : Check whether the input signal is long enough
% TO DO: determine length of the longest equivalent filter
% Do non-expansve transform if ext='per'
if(strcmp(flags.ext,'per'))
    doNoExt = 1;
else
    doNoExt = 0;
end


%% ----- step 3 : Run computation
treePath = nodesBForder(wt);
inLengths = nodeInLen(treePath,Ls,doNoExt,wt);
c = comp_wpfbt(f,wt.nodes(treePath),inLengths,'dec',flags.ext);