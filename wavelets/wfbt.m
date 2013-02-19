function c=wfbt(f,wt,varargin)
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
%   obtained from the |wfbtinit|_ function and modified arbitrarily or it
%   can be cell-array, which is used as a parameter in the internal call of
%   the |wfbtinit|_ function.
%
%   The following flag groups are supported:
%
%         'per','zpd','sym','symw','asym','asymw','ppd','sp0'
%                Type of the boundary handling.
%
%         'dwt','full'
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
%   A simple example of calling the |wfbt|_ function using the "full decomposition" wavelet tree:::
% 
%     f = gspi;
%     J = 7;
%     c = wfbt(f,{{'db',10},J},'full');
%     plotfwt(c,44100,90);
%
%   See also: iwfbt, wfbtinit
%


if(nargin<2)
   error('%s: Too few input parameters.',upper(mfilename));  
end

definput.import = {'fwt','wfbtcommon'};
%definput.flags.pack={'pack_null','pack'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Initialize the wavelet tree structure
wt = wfbtinit(wt,flags.treetype,'ana');

    
%% ----- step 1 : Verify f and determine its length -------
[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,upper(mfilename),0);
if(Ls<2)
   error('%s: Input signal seems not to be a vector of length > 1.',upper(mfilename));  
end

% Do non-expansve transform if ext='per'
if(strcmp(flags.ext,'per'))
    doNoExt = 1;
else
    doNoExt = 0;
end


%% ----- step 3 : Run computation
treePath = nodesBForder(wt);
inLengths = nodeInLen(treePath,Ls,doNoExt,wt);
rangeLoc = rangeInLocalOutputs(treePath,wt);
rangeOut = rangeInOutputs(treePath,wt);
c = comp_wfbt(f,wt.nodes(treePath),inLengths,rangeLoc,rangeOut,'dec',flags.ext);