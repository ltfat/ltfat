function outsig=operatoradj(Op,insig);
%OPERATORADJ  Apply the adjoint of an operator
%   Usage: c=operatoradj(Op,f);
%
%   `c=operatoradj(Op,f)` applies the adjoint operator of the operator *Op*
%   to the input signal *f*.  The operator object *Op* must have been
%   created using |operatornew|.
%
%   If *f* is a matrix, the transform will be applied along the columns
%   of *f*. If *f* is an N-D array, the transform will be applied along
%   the first non-singleton dimension.
%
%   See also: operatornew, operator, ioperator
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(Op)
  error('%s: First argument must be a operator definition structure.',upper(mfilename));
end;

switch(Op.type)
  case 'framemul'
    outsig=framemuladj(insig,Op.Fa,Op.Fs,Op.s);
  case 'spread'
    outsig=spreadadj(insig,Op.s);
end;

  
