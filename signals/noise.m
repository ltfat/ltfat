function outsig = noise(siglen,varargin)
% NOISE  Stochastic noise generator.
%   Usage: outsig = noise(siglen,nsigs,type);
%
%   Input parameters:
%       siglen    - Length of the noise (samples)
%       nsigs     - Number of signals (default is 1)
%       type      - type of noise: 'white', 'brown', 'pink', 'red'
%
%   Output parameters:
%       outsig      - siglen x nsigs signal vector
%
%   NOISE(siglen,nsigs) generates nsigs channels containing white noise of the
%   given type with the length of siglen. The signals are arranged as columns
%   in the output. If only siglen is given, a column vector is returned.
%
%   NOISE takes the following optional parameters:
%
%-     'white'  - Generate white (gaussian) noise. This is the default.
%
%-     'pink'   - Generate pink noise.
%
%-     'brown'  - Generate brown noise.
%
%-     'red'    - This is the same as brown noise.     
%
%   By default, the noise is normalized to have a unit energy, but this can
%   be changed by passing a flag to NORMALIZE.
%
%   See also: normalize

%   AUTHOR: Hagen Wierstorf and Peter L. Soendergaard.



% ------ Checking of input parameter -------------------------------------

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(siglen) || ~isscalar(siglen) || siglen<=0
    error('%s: siglen has to be a positive scalar.',upper(mfilename));
end

definput.import={'normalize'};
definput.importdefaults={'2'};
definput.flags.noisetype={'white','pink','brown','red'};
definput.keyvals.nsigs=1;
[flags,kv,nsigs]  = ltfatarghelper({'nsigs'},definput,varargin);

if flags.do_white
  outsig=randn(siglen,nsigs);
end;

if flags.do_brown || flags.do_red
  outsig=cumsum(randn(siglen,nsigs));
end;

if flags.do_pink
  % --- Handle trivial condition

  if siglen==1
    outsig=ones(1,nsigs);
    return;
  end;

  % ------ Computation -----------------------------------------------------
  fmax = floor(siglen/2)-1;
  f = (2:(fmax+1)).';
  % 1/f amplitude factor
  a = 1./sqrt(f);
  % Random phase
  p = randn(fmax,nsigs) + i*randn(fmax,nsigs);
  sig = repmat(a,1,nsigs).*p;

  outsig = ifftreal([ones(1,nsigs); sig; 1/(fmax+2)*ones(1,nsigs)],siglen);

end;

outsig=normalize(outsig,flags.norm);
