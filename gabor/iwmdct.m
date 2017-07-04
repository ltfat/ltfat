function [f,g]=iwmdct(c,g,Ls)
%IWMDCT  Inverse MDCT
%   Usage:  f=iwmdct(c,g);
%           f=iwmdct(c,g,Ls);
%
%   Input parameters:
%         c     : M*N array of coefficients.
%         g     : Window function.
%         Ls    : Final length of function (optional)
%   Output parameters:
%         f     : Input data
%
%   `iwmdct(c,g)` computes an inverse windowed MDCT with window *g*. The
%   number of channels is deduced from the size of the coefficient array *c*.
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |wilwin| for more details.
%
%   `iwmdct(f,g,Ls)` does the same, but cuts or zero-extends the final
%   result to length *Ls*.
%
%   `[f,g]=iwmdct(...)` additionally outputs the window used in the
%   transform. This is usefull if the window was generated from a
%   description in a string or cell array.
%
%   See also:  wmdct, wilwin, dgt, wildual, wilorth
%
%   References: prbr86 prjobr87 ma92 bohl96-1 

%   AUTHOR: Peter L. Søndergaard
%   TESTING: TEST_WMDCT

complainif_argnonotinrange(nargin,2,3,mfilename);

wasrow=0;
if isnumeric(g)
  if size(g,2)>1
    if size(g,1)>1
      error('g must be a vector');
    else
      % g was a row vector.
      g=g(:);
      
      % If the input window is a row vector, and the dimension of c is
      % equal to two, the output signal will also
      % be a row vector.
      if ndims(c)==2
        wasrow=1;
      end;
    end;
  end;
end;

M=size(c,1);
N=size(c,2);
W=size(c,3);

a=M;
L=M*N;

assert_L(L,0,L,a,2*M,'IWMDCT');

g=wilwin(g,M,L,'IWMDCT');

f=comp_idwiltiii(c,g);

% Check if Ls was specified.
if nargin==3
  f=postpad(f,Ls);
else
  Ls=L;
end;

f=comp_sigreshape_post(f,Ls,wasrow,[0; W]);

