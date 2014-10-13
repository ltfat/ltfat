function framegram(F,c,varargin)
%FRAMEGRAM  Easy visualization of energy in frame space
%   Usage: framegram(F,c,...);
%
%   `framegram(F,c)` plots the energy of the coefficients computed from
%   the input signal *x* using the frame *F* for analysis. This is
%   just a shorthand for::
%
%     plotframe(F,abs(frana(F,c)).^2);
%
%   Any additional arguments given to `framegram` are passed onto
%   |plotframe|.
%
%   See also: plotframe
    
plotframe(F,abs(frana(F,c)).^2,varargin{:});