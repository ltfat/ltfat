function h=framemul(f,Fa,Fs,sym,varargin)
%FRAMEMUL  Frame multiplier
%   Usage:  h=framemul(f,Fa,Fs,sym);
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
%   `framemul(f,Fa,Fs,sym)` applies the frame multiplier with symbol *sym*
%   to the signal *f*. The frame *Fa* is used for analysis and the frame
%   *Fs* for synthesis.
%
%   See also: iframemul, framemuladj
  
% Author: Peter L. Søndergaard

if nargin < 4
    error('%s: Too few input parameters.',upper(mfilename));
end;

h=frsyn(Fs,sym.*frana(Fa,f));






