function outsig = pinknoise(siglen,nsigs)
% PINKNOISE Generates a pink noise signal
%   Usage: outsig = pinknoise(siglen,nsigs);
%
%   Input parameters:
%       siglen    - Length of the noise (samples)
%       nsigs     - Number of signals (default is 1)
%
%   Output parameters:
%       outsig      - siglen x nsigs signal vector
%
%   PINKNOISE(siglen,nsigs) generates nsigs channels containing pink noise
%   (1/f spectrum) with the length of siglen. The signals are arranged as
%   columns in the output.

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input parameter -------------------------------------

error(nargchk(1,2,nargin));

if ~isnumeric(siglen) || ~isscalar(siglen) || siglen<=0
    error('%s: siglen has to be a positive scalar.',upper(mfilename));
end

if nargin==1
  nsigs=1;
end;

if ~isnumeric(nsigs) || ~isscalar(nsigs) || nsigs<=0
    error('%s: siglen has to be a positive scalar.',upper(mfilename));
end

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

% Scale output
%outsig = outsig ./ (max(abs(outsig(:)))+eps);
