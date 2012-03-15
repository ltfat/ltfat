function [tc,relres,iter,xrec] = framegrouplasso(F,x,lambda,varargin)
%FRAMEGROUPLASSO  Group LASSO regression in the TF-domain
%   Usage: [tc,xrec] = framegrouplasso(F,x,group,lambda,C,maxit,tol)
%
%   Input parameters:
%       F        : Frame definition
%       x        : Input signal
%       lambda   : Regularization parameter, controls sparsity of the
%                  solution
%   Output parameters:
%      tc        : Thresholded coefficients
%      relres    : Vector of residuals.
%      iter      : Number of iterations done.
%      xrec      : Reconstructed signal
%
%   `framegrouplasso(F,x)` solves the group LASSO regression problem in
%   the time-frequency domain: minimize a functional of the synthesis coefficients
%   defined as the sum of half the $l^2$ norm of the approximation error and
%   the mixed $l^1$ / $l^2$ norm of the coefficient sequence, with a penalization
%   coefficient lambda.
%  
%   The matrix of time-frequency coefficients is labelled in terms of groups and
%   members.  The obtained expansion is sparse in terms of groups, no
%   sparsity being imposed to the members of a given group. This is achieved
%   by a regularization term composed of $l^2$ norm within a group, and
%   $l^1$ norm with respect to groups.
%
%   `[tc,relres,iter] = framegrouplasso(...)` returns the residuals *relres* in
%   a vector and the number of iteration steps done, *maxit*.
%
%   `[tc,relres,iter,xrec] = framegrouplasso(...)` returns the reconstructed
%   signal from the coefficients, *xrec*. Note that this requires additional
%   computations.
%
%   The function takes the following optional parameters at the end of
%   the line of input arguments:
%
%     'freq'     Group in frequency (search for tonal components). This is the
%                default.
%
%     'time'     Group in time (search for transient components). 
%
%     'C',cval   Landweber iteration parameter: must be larger than
%                square of upper frame bound. Default value is the upper
%                frame bound.
%
%     'maxit',maxit
%                Stopping criterion: maximal number of iterations. Default value is 100.
%
%     'tol',tol  Stopping criterion: minimum relative difference between
%              norms in two consecutive iterations. Default value is
%              1e-2.
%
%     'print'    Display the progress.
%
%     'quiet'    Don't print anything, this is the default.
%
%     'printstep',p
%                If 'print' is specified, then print every p'th
%                iteration. Default value is 10;
%
%   The parameters *C*, *maxit* and *tol* may also be specified on the
%   command line in that order: `framegrouplasso(x,g,a,M,lambda,C,tol,maxit)`.
%
%   The solution is obtained via an iterative procedure, called Landweber
%   iteration, involving iterative group thresholdings.
%
%   The relationship between the output coefficients is given by ::
%
%     xrec = frsyn(F,tc);
%
%   See also: framelasso, framebounds

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isvector(x)
    error('Input signal must be a vector.');
end

% Define initial value for flags and key/value pairs.
definput.flags.group={'freq','time'};

definput.keyvals.C=[];
definput.keyvals.maxit=100;
definput.keyvals.tol=1e-2;
definput.keyvals.printstep=10;
definput.flags.print={'print','quiet'};

[flags,kv]=ltfatarghelper({'C','tol','maxit'},definput,varargin);

L=framelengthsignal(F,length(x));

if isempty(kv.C)
  [A_dummy,kv.C] = framebounds(F,L);
end;

tchoice=flags.do_time;
N = floor(length(x)/a);

% Normalization to turn lambda to a value comparable to lasso
if tchoice
    lambda = lambda * sqrt(N);
else
    lambda = lambda * sqrt(M);
end

% Various parameter initializations
threshold = lambda/kv.C;

% Initialization of thresholded coefficients
c0 = franaadj(F,x);
tc0 = c0;
relres = 1e16;
iter = 0;

% Main loop
while ((iter < kv.maxit)&&(relres >= kv.tol))
    tc = c0 - franaadj(F,frsyn(F,tc0));
    tc = tc0 + tc/kv.C;
    
    %  ------------ Convert to TF-plane ---------
    tc = framecoef2tf(F,tc);
    if tchoice
      tc = tc.';
    end;
    tc = groupthresh(tc,threshold,'soft');
    if tchoice
      tc=tc.';
    end;
    
    % Convert back from TF-plane
    tc=frametf2coef(F,tc);
    % -------------------------------------------
    
    relres = norm(tc(:)-tc0(:))/norm(tc0(:));
    tc0 = tc;
    iter = iter + 1;
    if flags.do_print
      if mod(iter,kv.printstep)==0        
        fprintf('Iteration %d: relative error = %f\n',iter,relres);
      end;
    end;
end

% Reconstruction
if nargout>3
  xrec = frsyn(F,tc);
end;

