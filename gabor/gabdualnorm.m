function [o1,o2]=gabdualnorm(gamma,g,a,M,L);
%GABDUALNORM  Measure of how close a window is to be a dual window. 
%   Usage:  dn=gabdualnorm(g,gamma,a,M);
%           dn=gabdualnorm(g,gamma,a,M,L);
%           [scal,res]=gabdualnorm(g,gamma,a,M);
%           [scal,res]=gabdualnorm(g,gamma,a,M,L);
%
%   Input parameters:
%         gamma  : input window..
%         g      : window function.
%         a      : Length of time shift.
%         M      : Number of modulations.
%         L      : Length of transform to consider
%   Output parameters:
%         dn     : dual norm.
%         scal   : Scaling factor
%         res    : Residual
%
%   GABDUALNORM(g,gamma,a,M) calculates how close gamma is to be a dual
%   window of the Gabor frame with window g and parameters _a and M.
%
%   The windows g and gamma may be vectors of numerical values, text strings
%   or cell arrays. See the help of GABWIN for more details.
%
%   GABDUALNORM(g,gamma,a,M,L) does the same, but considers a transform
%   length of L.
%
%   [scal,res]=GABDUALNORM(g,gamma,a,M) or
%   [scal,res]=GABDUALNORM(g,gamma,a,M,L) will compute two entities:
%   scal determines if the windows are scaled correctly, it must be 1 
%   for the windows to be dual. res is close to zero if the windows
%   (scaled correctly) are dual windows.
%
%   GABDUALNORM can be used to get the maximum relative reconstruction
%   error when using the two specified windows. Consider the following code
%   for some signal f, windows g, gamma, parameters _a and M and 
%   transform-length L (See help on DGT for how to obtain L):
%
%C     fr=idgt(dgt(f,g,a,M),gamma,a); 
%C     er=norm(f-fr)/norm(f);
%C     eest=gabdualnorm(g,gamma,a,M,L);
%
%   Then  _er < _eest for all possible input signals f.
%
%   To get a similar estimate for a tight window gt, simply use
%  
%C     eest=gabdualnorm(gt,gt,a,M,L);
%
%   See also:  gabframebounds, dgt

error(nargchk(4,5,nargin));

if nargin<5
  L=[];
end;

assert_squarelat(a,M,1,'GABDUALNORM',0);

if isnumeric(g)
  if ~isvector(g)
    error('%s: g must be a vector',upper(callfun));
  end;
  Ls=length(g);
else
  Ls=0;
end;

if isnumeric(gamma)
  if ~isvector(gamma)
    error('%s: gamma must be a vector',upper(callfun));
  end;
  Lwindow=length(gamma);
else
  Lwindow=0;
end;

[b,N,L]=assert_L(Ls,Lwindow,L,a,M,'GABDUALNORM');

[g,info_g]         = gabwin(g,a,M,L,'GABDUALNORM');
[gamma,info_gamma] = gabwin(gamma,a,M,L,'GABDUALNORM');
 
g=fir2long(g,L);
gamma=fir2long(gamma,L);

% Handle the Riesz basis (dual lattice) case.
if a>M

  % Calculate the right-hand side of the Wexler-Raz equations.
  rhs=dgt(gamma,g,a,M);
  scalconst=1;
  
else
  
  % Calculate the right-hand side of the Wexler-Raz equations.
  rhs=dgt(gamma,g,M,a);
  
  scalconst=a/M;
  
end;

if nargout<2
  % Subtract from the first element to make it zero, if the windows are
  % dual.
  rhs(1)=rhs(1)-scalconst;

  o1=norm(rhs(:),1);
else
  % Scale the first element to make it one, if the windows are dual.
  o1=rhs(1)/scalconst;
  o2=norm(rhs(2:end),1);
end;
  
