function [gamma]=wildual(g,M,L)
%WILDUAL  Canonical dual window.
%   Usage:  gamma=wildual(g,M);
%           gamma=wildual(g,M,L);
%
%   Input parameters:
%         g     : Gabor window.
%         M     : Number of modulations.
%         L     : Length of window. (optional)
%   Output parameters:
%         gamma : Canonical dual window.
%
%   WILDUAL(g,M) returns the dual window of the Wilson or WMDCT basis with
%   window g, parameter M and length equal to the length of the window g.
%
%   The window g may be a vector of numerical values, a text string or a
%   cell array. See the help of WILWIN for more details.
%
%   If the length of g is equal to 2*M, then the input window is assumed to
%   be a FIR window. In this case, the dual window also has length of
%   2*M. Otherwise the smallest possible transform length is chosen as the
%   window length.
%
%   WILDUAL(g,M,L) does the same, but now L is used as the length of the
%   Wilson basis. g will be cut or zero-extended to length L.
%
%   The input window g must be real and whole-point even. If g is not
%   whole-point even, then reconstruction using the dual window will not be
%   perfect. For a random window g, the window closest to g that satisfies
%   these restrictions can be found by
%
%C     g_wpe = real(g+involute(g));
%
%   All Gabor windows in the toolbox satisfies these restrictions unless
%   clearly stated otherwise.
%
%   See also:  dwilt, wilwin, wmdct, wilorth, isevenfunction

%   AUTHOR : Peter Soendergaard.
%   TESTING: TEST_DWILT
%   REFERENCE: OK

error(nargchk(2,3,nargin));
  
wasrow=0;

if size(g,2)>1
  if size(g,1)>1
    error('g must be a vector');
  else
    % g was a row vector.
    wasrow=1;
    g=g(:);
  end;
end;
  
assert_squarelat(M,M,1,'WILDUAL',0);

Lwindow=length(g);
Ls=Lwindow;

% There is no reason to do special support for FIR Wilson windows

if nargin<3
  [b,N,L]=assert_L(Ls,Lwindow,[],M,2*M,'WILDUAL');
else
  [b,N,L]=assert_L(L,Lwindow,L,M,2*M,'WILDUAL');
  g=fir2long(g,L);
end;

  
% If input is real, output must be real as well.
inputwasreal = isreal(g);

a=M;

gamma=2*comp_gabdual_long(g,a,2*M);

if inputwasreal
  gamma=real(gamma);
end;
      
if wasrow
  gamma=gamma.';
end;


