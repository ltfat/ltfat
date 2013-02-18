function h=framemuladj(f,Fa,Fs,sym,varargin)
%FRAMEMULADJ  Adjoint operator of frame multiplier
%   Usage: h=framemuladj(f,Fa,Fs,sym);
%
%   Input parameters:
%          Fa   : Analysis frame
%          Fs   : Synthesis frame
%          sym  : Symbol
%          f    : Input signal
%
%   Output parameters: 
%          h    : Output signal
%
%   `framemuladj(f,Fa,Fs,sym)` applies the adjoint of the frame multiplier
%   with symbol *sym* to the signal *f*. The frame *Fa* is used for analysis
%   and the frame *Fs* for synthesis.
%
%   `framemuladj(f,Fa,Fs,sym)` does the same using the frames *Fa* for
%   analysis and *Fs* for synthesis.
%
%   See also: framemul, iframemul
  
% Author: Peter L. Søndergaard

if nargin < 4
    error('%s: Too few input parameters.',upper(mfilename));
end;

% Swap the analysis and synthesis frames and conjugate the symbol.
h=frsyn(Fa,conj(sym).*frana(Fs,f));






