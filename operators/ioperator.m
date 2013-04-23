function outsig=ioperator(Op,insig);
%IOPERATOR  Apply inverse of operator
%   Usage: c=ioperator(Op,f);
%
%   `c=ioperator(Op,f)` applies the inverse op the operator *Op* to the
%   input signal *f*.  The operator object *Op* must have been created using
%   |operatornew|.
%
%   If *f* is a matrix, the transform will be applied along the columns
%   of *f*. If *f* is an N-D array, the transform will be applied along
%   the first non-singleton dimension.
%
%   See also: operatornew, operator, operatoradj
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(Op)
  error('%s: First agument must be a operator definition structure.',upper(mfilename));
end;

switch(Op.type)
  case 'framemul'
    outsig=iframemul(insig,Op.Fa,Op.Fs,Op.s);
  case 'spread'
    outsig=spreadinv(insig,Op.s);
end;

  
