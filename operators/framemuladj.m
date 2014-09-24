function h=framemuladj(f,Fa,Fs,s,varargin)
%FRAMEMULADJ  Adjoint operator of frame multiplier
%   Usage: h=framemuladj(f,Fa,Fs,s);
%
%   Input parameters:
%          Fa   : Analysis frame
%          Fs   : Synthesis frame
%          s    : Symbol
%          f    : Input signal
%
%   Output parameters: 
%          h    : Output signal
%
%   `framemuladj(f,Fa,Fs,s)` applies the adjoint of the frame multiplier
%   with symbol *s* to the signal *f*. The frame *Fa* is used for analysis
%   and the frame *Fs* for synthesis. This is equivalent to calling
%   `framemul(f,Fs,Fa,conj(s))`.
%
%
%   See also: framemul, iframemul
  
% Author: Peter L. Søndergaard

if nargin < 4
    error('%s: Too few input parameters.',upper(mfilename));
end;

% Swap the analysis and synthesis frames and conjugate the symbol.
h=frsyn(Fa,conj(s).*frana(Fs,f));






