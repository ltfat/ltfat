function outsig = pinknoise(varargin)
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
%
%   PINKNOISE is just a wrapper around NOISE(...,'pink');
%
%   See also: noise

outsig = noise(varargin{:},'pink');
