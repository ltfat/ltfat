function [g,a] = fwt2filterbank( w, J)
%FWT2FILTERBANK  FWT equivalent non-iterated filterbank
%   Usage: [g,a] = fwt2filterbank(wtdef)
%
%   Input parameters:
%          w : Wavelet filters definition
%          J : Number of filterbank iterations.
%
%   Output parameters:
%         g   : Cell array containing filters
%         a   : Vector of sub-/upsampling factors
%
%   `fwt2filterbank( w, J)` calculates the impulse responses *g* and the 
%   subsampling factors *a* of non-iterated filterbank, which is equivalent
%   to the wavelet filterbank tree described by *w* and *J*. The returned 
%   parameters can be used directly in |filterbank|, |ufilterbank|.
%   
%   The filters are scaled if *a* is not returned. 
%
%   The function is wrapper for calling |wfbt2filterbank|.
%   
%   See also: wfbtinit, wfbt2filterbank, filterbank


if(nargin<2)
    error('%s: Not enough input arguments',upper(mfilename));
end

if nargout<2
  g = wfbt2filterbank({w,J,'dwt'});
elseif nargout == 2
  [g,a] = wfbt2filterbank({w,J,'dwt'});   
end

