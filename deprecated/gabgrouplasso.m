function [tc,relres,iter,xrec] = gabgrouplasso(x,g,a,M,lambda,varargin)
%GABGROUPLASSO  Group LASSO regression in Gabor domain
%   Usage: [tc,xrec] = gabgrouplasso(x,g,a,M,group,lambda,C,maxit,tol)
%
%   `gabgrouplasso` has been deprecated. Please use |framegrouplasso|_ instead.
%
%   A call to `gabgrouplasso(x,g,a,M,lambda)` can be replaced by ::
%
%     F=newframe('dgt',[],g,a,M);
%     tc=framegrouplasso(F,lambda);
%
%   Any additional parameters passed to `gabgrouplasso` can be passed to
%   |framegrouplasso|_ in the same manner.
%
%   See also: newframe, framegrouplasso

warning(['LTFAT: GABGROUPLASSO has been deprecated, please use FRAMEGROUPLASSO ' ...
         'instead. See the help on GABGROUPLASSO for more details.']);   

F=newframe('dgt',[],g,a,M);
[tc,relres,iter,xrec] = framegrouplasso(F,lambda,varargin{:});
