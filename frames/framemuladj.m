function h=framemuladj(f,varargin)
%FRAMEMULADJ  Adjoint operator of frame multiplier
%   Usage: h=framemuladj(f,F,sym);
%         h=framemuladj(f,Fa,Fs,sym);
%
%   Input parameters:
%          F    : Frame used for both analysis and synthesis
%          Fa   : Analysis frame
%          Fs   : Synthesis frame
%          sym  : Symbol
%          f    : Input signal
%
%   Output parameters: 
%          h    : Output signal
%
%   `framemuladj(f,F,sym)` applies the adjoint of the frame multiplier with
%   symbol *sym* to the signal *f*. The frame *F* is used for both analysis
%   and synthesis.
%
%   `framemul(f,Fa,Fs,sym)` does the same using the frames *Fa* for
%   analysis and *Fs* for synthesis.
%
%   See also: framemul
  
% Author: Peter L. Søndergaard

if nargin < 2
    error('%s: Too few input parameters.',upper(mfilename));
end;

if nargin==3
  Fa=varargin{1};
  Fs=varargin{1};
  sym=varargin{2};
end;

if nargin==4
  Fa=varargin{1};
  Fs=varargin{2};
  sym=varargin{3};
end;

% Swap the analysis and synthesis frames and conjugate the symbol.
h=frsyn(Fa,conj(sym).*frana(Fs,f));






