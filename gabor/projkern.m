function c=projkern(c,p2,p3,p4,p5);
%PROJKERN  Projection onto generating kernel space
%   Usage:  cout=projkern(cin,a);
%           cout=projkern(cin,g,a);
%           cout=projkern(cin,ga,gs,a);
%
%   Input parameters:
%         cin   : Input coefficients
%         g     : analysis/synthesis window
%         ga    : analysis window
%         gs    : synthesis window
%         a     : Length of time shift.
%   Output parameters:
%         cout  : Output coefficients
%
%   cout=PROJKERN(cin,a) projects the symbol c of a Gabor multiplier
%   onto the space of realisiable symbols. A tight window generated from a
%   Gaussian will be used.
%
%   The rationale for this function is a follows: Because the coefficient
%   space of a Gabor frame is larger than the signal space (since the frame
%   is redundant) then there are many coefficients that correspond to the
%   same signal. Therefore, when designing the symbol of a gabor multiplier,
%   the choice of symbol is more limited than one would expect.
%
%   Therefore, you might desire to work with the symbol cin, but you in
%   are really working with cout.
%
%   cout=PROJKERN(cin,g,a) does the same, using the window g for analysis
%   and synthesis.
%
%   cout=PROJKERN(cin,ga,gs,a) does the same, but for different analysis
%   ga and synthesis gs windows.
%
%   See also: gabmul
%
%   Demos: demo_gabmul

error(nargchk(2,4,nargin));

M=size(c,1);
N=size(c,2);

if nargin==2
  a=p2;
  L=a*N;
  ga=gabtight(a,M,L);
  gs=ga;
end;

if nargin==3;
  ga=p2;
  gs=p2;
  a=p3;
  L=a*N;
end;

if nargin==4;  
  ga=p2;
  gs=p3;
  a=p4;
  L=a*N;
end;

assert_squarelat(a,M,1,'PROJKERN');

c=dgt(idgt(c,gs,a),ga,a,M);

