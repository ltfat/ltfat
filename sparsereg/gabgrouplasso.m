function [tc,xrec,rel_diff,iter] = gabgrouplasso(x,g,a,M,lambda,varargin)
%GABGROUPLASSO  Group LASSO regression in Gabor domain
%   Usage: [tc,xrec] = gabgrouplasso(x,g,a,M,group,lambda,C,itermax,tol)
%   Input parameters:
%       x        : Input signal
%       g        : Synthesis window function
%       a        : Length of time shift
%       M        : Number of channels
%       lambda   : Regularization parameter, controls sparsity of the
%                   solution
%       group    : 'time' or 'freq' (default)
%       C        : Landweber iteration parameter: must be larger than
%                   square of upper frame bound.
%       itermax  : Stopping criterion: maximal number of iterations
%       tol      : Stopping criterion: minimum relative difference between
%                   norms in two consecutive iterations
%   Output parameters:
%      tc        : Thresholded coefficients
%      xrec      : Reconstructed signal
%
%   GABGROUPLASSO(x,g,a,M) solves the group LASSO
%   regression problem in the Gabor domain: minimize a functional
%   of the synthesis coefficients defined as the sum of half the l2 norm of 
%   the approximation error and the mixed l1/l2 norm of the coefficient
%   sequence, with a penalization coefficient lambda.
%  
%   The matrix of Gabor coefficients is labelled in terms of groups and
%   members.  The obtained expansion is sparse in terms of groups, no
%   sparsity being imposed to the members of a given group. This is achieved
%   by a regularization term composed of l2 norm within a group, and l1 norm
%   with respect to groups.
%
%   The function takes the following optional parameters at the end of
%   the line of input arguments:
%
%-      'freq' - Group in frequency (search for tonal components). This is the
%                default.
%
%-      'time' - Group in time (search for transient components). 
%
%-      'C',cval - Landweber iteration parameter: must be larger than
%                square of upper frame bound. Default value is the upper
%                frame bound.
%
%-      'itermax',itermax - Stopping criterion: maximal number of
%                iterations. Default value is 100.
%
%       'tol',tol - Stopping criterion: minimum relative difference between
%                norms in two consecutive iterations. Default value is
%                1e-2.
%
%   The parameters C, itermax and tol may also be specified on the
%   command line in that order: GABGROUPLASSO(x,g,a,M,lambda,C,itermax,tol).
%
%   The solution is obtained via an iterative procedure, called Landweber
%   iteration, involving iterative group thresholdings.
%
%   The relationship between the output coefficients is given by
%
%C      xrec = idgt(tc,g,a);
%
%   See also: gablasso, gabframebounds

if nargin<5
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isvector(x)
    error('Input signal must be a vector.');
end

% Define initial value for flags and key/value pairs.
defnopos.flags.group={'freq','time'};

defnopos.keyvals.C=[];
defnopos.keyvals.itermax=100;
defnopos.keyvals.tol=1e-2;

[flags,kv]=ltfatarghelper({'C','itermax','tol'},defnopos,varargin);

% Determine transform length, and calculate the window.
[x,g,L] = gabpars_from_windowsignal(x,g,a,M,[],'GABGROUPLASSO');

if isempty(kv.C)
  [A_dummy,kv.C] = gabframebounds(g,a,M,L);
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
c0 = dgt(x,g,a,M);
tc0 = c0;
rel_diff = 1e16;
iter = 0;

% Main loop
while ((iter < kv.itermax)&&(rel_diff >= kv.tol))
    tc = c0 - dgt(idgt(tc0,g,a),g,a,M);
    tc = tc0 + tc/kv.C;
    if tchoice
        tc = tc';
    end;
    tc = groupthresh(tc,threshold,'soft');
    if tchoice
        tc=tc';
    end;
    rel_diff = norm(tc(:)-tc0(:))/norm(tc0(:));
    tc0 = tc;
    iter = iter + 1;
    fprintf('Iteration %d: relative error = %f\n',iter,rel_diff);
end

% Reconstruction
if nargout>1
  xrec = idgt(tc,g,a);
end;
