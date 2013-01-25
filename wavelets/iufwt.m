function f = iufwt(c,g,J,varargin)
%IUFWT   Inverse Undecimated Fast Wavelet Transform 
%   Usage:  f = iufwt(c,g,J,...)
%
%   Input parameters:
%         c     : Coefficients stored in a cell-array.
%         g     : Synthesis wavelet filters.
%         J     : Number of filterbank iterations.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   `f = iufwt(c,g,J)` reconstruct signal *f* from the wavelet
%   coefficients *c* using the wavelet filterbank consisting of the *J* levels
%   of the basic synthesis filterbank defined by *g* using the "a-trous" algorithm.
%
%   The following flags are supported:
%
%         'per','zpd','sym','symw','asym','asymw','ppd','sp0'
%                Type of the boundary handling.
%
%   Node that the same flag as in the `ufwt` function have to be used.
%
%   Please see the help on |ufwt|_ for a description of the parameters.
%
%   See also:  ufwt, fwtinit
%
%   References: ma98


if nargin<3
   error('%s: Too few input parameters.',upper(mfilename));
end;

% Initialize the wavelet filters structure
g = fwtinit(g,'syn');

%% PARSE INPUT
definput.keyvals.Ls=[];    
definput.import = {'fwt'};

if(iscell(c))
    [flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);
else
    error('%s: Unrecognized coefficient format.',upper(mfilename));
end

% Estimate output signal length from the number of coefficients
[sigHalfLen,W] = size(c{end});
if(strcmp(flags.ext,'per'))
   Ls = sigHalfLen;  
else
   Ls = sigHalfLen - (length(g.g{end})-1);
end

f = comp_ifwt_all(c,g.g,J,g.a,Ls,'undec',flags.ext);







