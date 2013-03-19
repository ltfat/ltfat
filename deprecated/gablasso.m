function [tc,relres,iter,xrec] = gablasso(x,g,a,M,lambda,varargin)
%GABLASSO  LASSO regression in Gabor domain
%   Usage: [tc,xrec] = gablasso(x,a,M,lambda,C,tol,maxit)
%
%   `gablasso` has been deprecated. Please use |framelasso| instead.
%
%   A call to `gablasso(x,g,a,M,lambda)` can be replaced by ::
%
%     F=newframe('dgt',[],g,a,M);
%     tc=framelasso(F,lambda);
%
%   Any additional parameters passed to `gablasso` can be passed to
%   |framelasso| in the same manner.
%
%   See also: newframe, framelasso

warning(['LTFAT: GABLASSO has been deprecated, please use FRAMELASSO ' ...
         'instead. See the help on GABLASSO for more details.']);   

F=newframe('dgt',[],g,a,M);
[tc,relres,iter,xrec] = framelasso(F,lambda,varargin{:});
