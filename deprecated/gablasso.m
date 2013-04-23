function [tc,relres,iter,xrec] = gablasso(x,g,a,M,lambda,varargin)
%GABLASSO  LASSO regression in Gabor domain
%   Usage: [tc,xrec] = gablasso(x,a,M,lambda,C,tol,maxit)
%
%   `gablasso` has been deprecated. Please use |franalasso| instead.
%
%   A call to `gablasso(x,g,a,M,lambda)` can be replaced by ::
%
%     F=frame('dgt',[],g,a,M);
%     tc=franalasso(F,lambda);
%
%   Any additional parameters passed to `gablasso` can be passed to
%   |franalasso| in the same manner.
%
%   See also: frame, franalasso

warning(['LTFAT: GABLASSO has been deprecated, please use FRANALASSO ' ...
         'instead. See the help on GABLASSO for more details.']);   

F=newframe('dgt',[],g,a,M);
[tc,relres,iter,xrec] = franalasso(F,lambda,varargin{:});
