function V=latticetype2matrix(L,a,M,lt);
%LATTICETYPE2MATRIX  Convert lattice description to matrix form
%   Usage: V=latticetype2matrix(L,a,M,lt);
%
%   `V=latticetype2matrix(L,a,M,lt)` converts a standard description of a
%   lattice using the *a*, *M* and *lt* parameters into a $2\times 2$
%   integer matrix description. The conversion is *only* valid for the
%   specified transform length *L*.
%
%   The output will be in lower triangular Hemite normal form.
%
%   For more information, see
%   `<http://en.wikipedia.org/wiki/Hermite_normal_form>`_.
%
%   An example:::
%
%     V = latticetype2matrix(120,10,12,[1 2])
%
%   See also: matrix2latticetype

L2=dgtlength(L,a,M,lt);

if L~=L2
    error('%s: Invalid transform length.',upper(mfilename));
end;

b=L/M;
s=b/lt(2)*lt(1);
V=[a 0;...
   s b];

