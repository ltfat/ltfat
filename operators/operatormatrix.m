function T=operatormatrix(Op)
%OPERATORMATRIX  Matrix representation of an operator
%   Usage: T=operatormatrix(Op);
%
%   `T=operatormatrix(Op)` returns the matrix representation *T* of the
%   operator *Op*. The operator object *Op* must have been created using
%   |operatornew|.
%
%   See also: operatornew, operator, operatoreigs

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(Op)
  error('%s: First argument must be a operator definition structure.',upper(mfilename));
end;

T=operator(Op,eye(Op.L));
