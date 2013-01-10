function [V,D]=framemuleigs(Fa,Fs,coef,varargin)
%FRAMEMULEIGS  Eigenpairs of frame multiplier
%   Usage:  h=framemuleigs(F,c,K);
%           h=framemuleigs(F,c,K,...);
%
%   Input parameters:
%         Fa    : Frame analysis definition
%         Fs    : Frame analysis definition
%         K     : Number of eigenvectors to compute.
%         c     : symbol of Gabor multiplier
%   Output parameters:
%         V     : Matrix containing eigenvectors.
%         D     : Eigenvalues.
%
%   `framemuleigs(F,c,K)` computes the *K* largest eigenvalues and eigen-
%   vectors of the frame multiplier with symbol *c*.
%
%   If *K* is empty, then all eigenvalues/pairs will be returned.
%
%   `framemuleigs` takes the following parameters at the end of the line of input
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
%   Examples:
%   ---------
%
%   The following example calculates and plots the first eigenvector of the
%   Gabor multiplier given by the |batmask|_ function. Note that the mask
%   must be converted to a column vector to work with in this framework:::
%
%     mask=batmask;
%     [Fa,Fs]=framepair('dgt','gauss','dual',10,40);
%     [V,D]=framemuleigs(Fa,Fs,mask(:));
%     sgram(V(:,1),'dynrange',90);
%
%   See also: gabmul, dgt, idgt, gabdual, gabtight

% Change this to 1 or 2 to see the iterative method in action.
printopts=0;

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if nargout==2
  doV=1;
else
  doV=0;
end;

definput.keyvals.K=6;
definput.keyvals.maxit=100;
definput.keyvals.tol=1e-9;
definput.keyvals.crossover=200;
definput.flags.print={'quiet','print'};
definput.flags.method={'auto','iter','full'};


[flags,kv,K]=ltfatarghelper({'K'},definput,varargin);

% Do the computation. For small problems a direct calculation is just as
% fast.

L=framelengthcoef(Fa,size(coef,1));

if (flags.do_iter) || (flags.do_auto && L>kv.crossover)
    
  if flags.do_print
    opts.disp=1;
  else
    opts.disp=0;
  end;
  opts.isreal = false;
  opts.maxit  = kv.maxit;
  opts.tol    = kv.tol;
  
  if doV
    [V,D] = eigs(@(x) frsyn(Fs,coef.*frana(Fa,x)),L,K,'LM',opts);
  else
    D     = eigs(@(x) frsyn(Fs,coef.*frana(Fa,x)),L,K,'LM',opts);
  end;

else
  % Compute the transform matrix.
  bigM=frsyn(F,diag(coef)*frana(F,eye(L)));

  if doV
    [V,D]=eig(bigM);
  else
    D=eig(bigM);
  end;

end;

% The output from eig and eigs is a diagonal matrix, so we must extract the
% diagonal.
D=diag(D);

% Sort them in descending order
[~,idx]=sort(abs(D),1,'descend');
D=D(idx(1:K));

if doV
  V=V(:,idx(1:K));
end;

% Clean the eigenvalues, if we know that they are real-valued
%if isreal(ga) && isreal(gs) && isreal(c)
%  D=real(D);
%end;
