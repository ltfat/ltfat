function outsig=operator(Op,insig);
%OPERATOR  Apply operator
%   Usage: c=operator(Op,f);
%
%   `c=operator(Op,f)` applies the operator *Op* to the input signal *f*.
%   The operator object *Op* must have been created using |operatornew|.
%
%   If *f* is a matrix, the transform will be applied along the columns
%   of *f*. If *f* is an N-D array, the transform will be applied along
%   the first non-singleton dimension.
%
%   See also: operatornew, ioperator, operatoradj
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(Op)
  error('%s: First agument must be a operator definition structure.',upper(mfilename));
end;

switch(Op.type)
  case 'framemul'
    outsig=framemul(insig,Op.Fa,Op.Fs,Op.s);
  case 'spread'
    outsig=spreadop(insig,Op.s);
end;

  
