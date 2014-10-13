function Opout=operatorappr(Op,T)
%OPERATORAPPR  Best approximation by operator
%   Usage: c=operatorappr(Op,K);
%
%   `Opout=operatorappr(Opin,T)` computes the an operator *Opout* of the
%   same type as *Opin* that best approximates the matrix *T* in the
%   Frobenious norm of the matrix (the Hilbert-Schmidt norm of the
%   operator).
%
%   For some operator classes, the approximation is always exact, so that
%   `operator(Opout,f)` computes the exact same result as `T'*f`.
%
%   See also: operatornew, operator, operatoreigs
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(Op)
  error('%s: First argument must be a operator definition structure.',upper(mfilename));
end;

switch(Op.type)
  case 'framemul'
    s=framemulappr(Op.Fa,Op.Fs,T);
    Opout=operatornew('framemul',Op.Fa,Op.Fs,s);
  case 'spread'
    s=spreadfun(T);
    Opout=operatornew('spread',s);
end;

  
