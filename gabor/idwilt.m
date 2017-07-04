function [f,g]=idwilt(c,g,Ls)
%IDWILT  Inverse discrete Wilson transform
%   Usage:  f=idwilt(c,g);
%           f=idwilt(c,g,Ls);
%
%   Input parameters:
%      c     : $2M \times N$ array of coefficients.
%      g     : Window function.
%      Ls    : Final length of function (optional)
%   Output parameters:
%      f     : Input data
%
%   `idwilt(c,g)` computes an inverse discrete Wilson transform with window *g*.
%   The number of channels is deduced from the size of the coefficient array *c*.
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |wilwin| for more details.
%  
%   `idwilt(f,g,Ls)` does the same, but cuts of zero-extend the final
%   result to length *Ls*.
%
%   `[f,g]=idwilt(...)` additionally outputs the window used in the
%   transform. This is usefull if the window was generated from a
%   description in a string or cell array.
%
%   See also:  dwilt, wilwin, dgt, wilorth
%
%   References: bofegrhl96-1 liva95

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: TEST_DWILT
%   REFERENCE: OK

complainif_argnonotinrange(nargin,2,3,mfilename);

M=size(c,1)/2;
N=2*size(c,2);
W=size(c,3);

a=M;
L=M*N;

assert_L(L,0,L,a,2*M,'IDWILT');

[g,info]=wilwin(g,M,L,'IDWILT');

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

