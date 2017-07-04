function coef=spreadfun(T)
%SPREADFUN  Spreading function of a matrix
%   Usage:  c=spreadfun(T);
%
%   `spreadfun(T)` computes the spreading function of the operator *T*,
%   represented as a matrix. The spreading function represent the operator *T*
%   as a weighted sum of time-frequency shifts. See the help text for
%   |spreadop| for the exact definition.
%
%   See also:  spreadop, tconv, spreadinv, spreadadj

complainif_argnonotinrange(nargin,1,1,mfilename);

if ndims(T)>2 || size(T,1)~=size(T,2)
    error('Input symbol T must be a square matrix.');
end;

L=size(T,1);

% The 'full' appearing on the next line is to guard the mex file.
coef=comp_col2diag(full(T));

coef=fft(coef)/L;

