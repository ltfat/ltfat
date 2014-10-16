function [tc,relres,iter,xrec] = gabgrouplasso(x,g,a,M,lambda,varargin)
%GABGROUPLASSO  Group LASSO regression in Gabor domain
%   Usage: [tc,xrec] = gabgrouplasso(x,g,a,M,group,lambda,C,maxit,tol)
%
%   `gabgrouplasso` has been deprecated. Please use |franagrouplasso| instead.
%
%   A call to `gabgrouplasso(x,g,a,M,lambda)` can be replaced by ::
%
%     F=frame('dgt',g,a,M);
%     tc=franagrouplasso(F,lambda);
%
%   Any additional parameters passed to `gabgrouplasso` can be passed to
%   |franagrouplasso| in the same manner.
%
%   See also: frame, franagrouplasso

warning(['LTFAT: GABGROUPLASSO has been deprecated, please use FRANAGROUPLASSO ' ...
          'instead. See the help on FRANAGROUPLASSO for more details.']);  

F=newframe('dgt',[],g,a,M);
[tc,relres,iter,xrec] = framegrouplasso(F,lambda,varargin{:});
