function out=spreadinv(p1,p2);
%SPREADINV  Apply inverse spreading operator
%   Usage: h=spreadinv(f,c);
%
%   `spreadinv(c)` computes the symbol of the inverse of the spreading
%   operator with symbol *c*.
%
%   `spreadinv(f,c)` applies the inverse of the spreading operator with
%   symbol *c* to the input signal *f*.
%
%   See also: spreadfun, tconv, spreadfun, spreadadj

complainif_argnonotinrange(nargin,1,2,mfilename);

% FIXME This function should handle sparse symbols, and use a iterative
% method instead of creating the full symbol.

% FIXME This function should handle f though comp_reshape_pre and post.
if nargin==1
  coef=p1;
else
  f=p1;
  coef=p2;
end;

if ndims(coef)>2 || size(coef,1)~=size(coef,2)
    error('Input symbol T must be a square matrix.');
end;

L=size(coef,1);

% Create a matrix representation of the operator.
coef=ifft(full(coef))*L;
T=comp_col2diag(coef);

if nargin==1
  
  % Calculate the inverse symbol.
  out=spreadfun(inv(T));
    
else
  
  out=T\f;
  
end;
