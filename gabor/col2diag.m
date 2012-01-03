function cout=col2diag(cin)
%COL2DIAG  Move columns of a matrix to diagonals
%   Usage:  cout=col2diag(cin);
%
%   COL2DIAG(cin) will rearrange the elements in the square matrix cin so
%   that columns of cin appears as diagonals. Column number _n will appear
%   as diagonal number _-n and _L-n, where _L is the size of the matrix.
%
%   The function is its own inverse.
%
%   COL2DIAG performs the underlying coordinate transform for spreading
%   function and Kohn-Nirenberg calculus in the finite, discrete setting.
%M
%   See also: spreadop, spreadfun, tconv
  
%   AUTHOR : Peter Soendergaard.
%   TESTING: TEST_SPREAD
%   REFERENCE: OK

% Assert correct input.
error(nargchk(1,1,nargin));

if ndims(cin)~=2 || size(cin,1)~=size(cin,2)
  error('Input matrix must be square.');
end;

if ~isnumeric(cin)
  error('Input must be numerical.');
end;

cout=comp_col2diag(full(cin));
%OLDFORMAT
