function V=latticetype2matrix(L,a,M,lt);
%LATTICETYPE2MATRIX  Convert lattice description to matrix form
%   Usage: V=latticetype2matrix(L,a,M,lt);
%
%   `V=latticetype2matrix(L,a,M,lt)` converts a standard description of a
%   lattice using the *a*, *M* and *lt* parameters into a $2\times 2$
%   integer matrix description. The conversion is *only* valid for the
%   specified transform length *L*.
%
%   The output will be in Hemite normal form, see |hermitenf|_.
%
%   An example:::
%
%     V = latticetype2matrix(120,10,12,[1 2])
%
%   See also: matrix2latticetype, hermitenf

L2=dgtlength(L,a,M,lt);

if L~=L2
    error('Invalid transform length.');
end;

b=L/M;
s=b/lt(2)*lt(1);
V=[a 0;...
   s b];

