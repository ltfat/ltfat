function f = ifwt(c,g,J,varargin)
%IFWT   Inverse Fast Wavelet Transform 
%   Usage:  f = ifwt(c,g,J)
%           f = ifwt(c,g,J,Ls,...)
%
%   Input parameters:
%         c     : Coefficients stored in J+1 cell-array.
%         g     : Synthesis wavelet filters.
%         J     : Number of filterbank iterations.
%         Ls    : Length of the reconstructed signal.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   `f = ifwt(c,g,J)` reconstructs signal *f* from the wavelet coefficients
%   *c* using *J*-iteration synthesis filter bank build from the basic synthesis
%   filterbank defined by *g*. The fast wavelet transform algorithm 
%   (or Mallat's algorithm) is employed. 
%
%   The supported boundary conditions are:
%
%     'per'    Periodic boundary extension. This is the default.
%
%     'zpd'    Zeros are considered outside of the signal (coefficient) support. 
%
%     'sym'    Half-point symmetric extension.
%
%     'symw'   Whole-point symmetric extension
%
%     'asym'   Half-point antisymmetric extension
%
%     'asymw'  Whole point antisymmetric extension
%
%     'ppd'    Periodic padding, same as `'per'` but the result is expansive representation
%
%     'sp0'    Repeating boundary sample
%
%   Node that the same flag as in the |fwt|_ function have to be used.
%
%   Please see the help on |fwt|_ for a description of the parameters.
%
%   Examples:
%   ---------
%   
%   A simple example showing perfect reconstruction:::
% 
%     f = gspi;
%     J = 8;
%     c = fwt(f,{'db',8},J);
%     fhat = ifwt(c,{'db',8},J,length(f));
%     % The following should give (almost) zero
%     norm(f-fhat)
%   
%   See also:  fwt, fwtinit
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


if isempty(Ls)
    % Estimate output signal length from the number of coefficients
   [sigHalfLen,W] = size(c{end});
       if(strcmp(flags.ext,'per'))
            % estimated Ls can be one sample more, if the original input
            % signal length was odd
           Ls = g.a(end)*sigHalfLen; 
       else
            % estimated Ls can be one sample more, if the original input
            % signal length plus length(h{1})-1 was an even number
           Ls = g.a(end)*sigHalfLen - (length(g.filts{end}.h)-2);
       end
end


f = comp_ifwt_all(c,g.filts,J,g.a,Ls,'dec',flags.ext);


