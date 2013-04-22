function outsig=operatoreigs(Op,K,varargin);
%OPERATOREIGS  Apply the adjoint of an operator
%   Usage: c=operatoreigs(Op,K);
%
%   `[V,D]=operatoreigs(Op,K)` computes the *K* largest eigenvalues and
%   eigenvectors of the operator *Op* to the input signal *f*.  The operator
%   object *Op* must have been created using |operatornew|.
%
%   If *K* is empty, then all eigenvalues/pairs will be returned.
%
%   `D=operatoreigs(...)` computes only the eigenvalues.
%
%   `operatoreigs` takes the following parameters at the end of the line of input
%   arguments:
%
%     'tol',t      Stop if relative residual error is less than the
%                  specified tolerance. Default is 1e-9 
%
%     'maxit',n    Do at most n iterations.
%
%     'iter'       Call `eigs` to use an iterative algorithm.
%
%     'full'       Call `eig` to solve the full problem.
%
%     'auto'       Use the full method for small problems and the
%                  iterative method for larger problems. This is the
%                  default. 
%
%     'crossover',c
%                  Set the problem size for which the 'auto' method
%                  switches. Default is 200.
%
%     'print'      Display the progress.
%
%     'quiet'      Don't print anything, this is the default.
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
    outsig=framemuleigs(Op.Fa,Op.Fs,Op.s,K,varargin{:});
  case 'spread'
    outsig=spreadeigs(K,Op.s);
end;

  
