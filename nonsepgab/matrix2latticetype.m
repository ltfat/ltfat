function [a,M,lt] = matrix2latticetype(L,V);
%MATRIX2LATTICETYPE  Convert matrix form to standard lattice description
%   Usage: [a,M,lt] = matrix2latticetype(L,V);
%
%   `[a,M,lt]=matrix2latticetype(L,V)` converts a $2\times 2$ integer matrix
%   description into the standard description of a lattice using the *a*,
%   *M* and *lt*. The conversion is *only* valid for the specified transform
%   length *L*.
%
%   An example:::
%
%     [a,M,lt] = matrix2latticetype(120,[10 0; 5 10])
%
%   See also: latticetype2matrix, hermitenf

V=hermitenf(V);

a=V(1,1);
b=V(2,2);
s=V(2,1);
M=L/b;

k=gcd(s,b);
lt=[s/k b/k];

