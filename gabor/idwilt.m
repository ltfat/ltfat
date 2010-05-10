function [f]=idwilt(c,g,Ls)
%IDWILT  Inverse discrete Wilson transform.
%   Usage:  f=idwilt(c,g);
%           f=idwilt(c,g,Ls);
%
%   Input parameters:
%         c     : M*N array of coefficients.
%         g     : Window function.
%         Ls    : Final length of function (optional)
%   Output parameters:
%         f     : Input data
%
%   IDWILT(c,g) computes an inverse discrete Wilson transform with window g.
%   The number of channels is deduced from the size of the coefficient array c.
%
%   The window g may be a vector of numerical values, a text string or a
%   cell array. See the help of WILWIN for more detailts.
%  
%   IDWILT(f,g,Ls) does the same, but cuts of zero-extend the final
%   result to length Ls.
%
%   See also:  dwilt, wilwin, dgt, wilorth
%
%R  bofegrhl96-1 liva95

%   AUTHOR : Peter Soendergaard.
%   TESTING: TEST_DWILT
%   REFERENCE: OK

error(nargchk(2,3,nargin));

M=size(c,1)/2;
N=2*size(c,2);
W=size(c,3);

a=M;
L=M*N;

assert_L(L,0,L,a,2*M,'IDWILT');

[g,info]=comp_window(g,a,2*M,L,1,'IDWILT');

wasrow=0;
if (ndims(c)==2 && info.wasrow)
  wasrow=1;
end;

f=comp_idwilt(c,g);

% Check if Ls was specified.
if nargin==3
  f=postpad(f,Ls);
else
  Ls=L;
end;

f=comp_sigreshape_post(f,Ls,wasrow,[0; W]);

    
