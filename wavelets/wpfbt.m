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
%   `c=wpfbt(f,wt)` returns coefficients *c* obtained applying wavelet filterbank tree
%   defined by *wt* to the input data *f*. If *f* is a matrix, the transformation 
%   is applied to each of *W* columns. 
%
%   
%
%
%
%
%   The following flag groups are supported:
%
%         'per','zpd','sym','symw','asym','asymw','ppd','sp0'
%                Type of the boundary handling.
%
%   Please see the help on |fwt|_ for a description of the flags.
%
%   Examples:
%   ---------
%   
%   A simple example of calling the |wfbt|_ function using the "full decomposition" wavelet tree:::
% 
%     f = gspi;
%     J = 6;
%     c = wfbt(f,{{'db',10},J},'full');
%     plotfwt(c);
%
%   See also: iwfbt, wfbtinit


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

%% ----- step 2 : Check whether the input signal is long enough
% TO DO: determine length of the longest equivalent filter

%% ----- step 3 : Run computation
c = comp_wfbt(f,wt,'dec',flags.ext);