function inoutsig = gaindb(inoutsig,gn)
%GAINDB  Increase/decrease level of signal
%   Usage:  outsig = gaindb(insig,gn);
%
%   GAINDB(insig,gn) changes the energy level of the signal by gn dB.
%
%   See also: rms

%   AUTHOR: Peter L. Soendergaard, 2009

% ------ Checking of input parameters ---------
  
error(nargchk(2,2,nargin));

if ~isnumeric(inoutsig)
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(gn) || ~isscalar(gn) 
  error('%s: g must be a scalar.',upper(mfilename));
end;

% ------ Computation --------------------------

inoutsig = inoutsig*10^(gn/20);
