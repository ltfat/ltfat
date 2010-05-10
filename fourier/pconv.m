function h=pconv(p1,p2,p3)
%PCONV  Periodic convolution
%   Usage:  h=pconv(f,g)
%           h=pconv(ptype,f,g); 
%
%   PCONV(f,g) computes the periodic convolution of f and g. The convolution
%   is given by
%
%M              L-1
%M     h(l+1) = sum f(k+1) * g(l-k+1)
%M              k=0
%
%F  \[h\left(l+1\right)=\sum_{k=0}^{L-1}f\left(k+1\right)g\left(l-k+1\right)\]
%   PCONV('r',f,g) computes the alternative where g is reversed given by
%
%M              L-1
%M     h(l+1) = sum f(k+1) * conj(g(k-l+1))
%M              k=0
%
%F  \[h\left(l+1\right)=\sum_{k=0}^{L-1}f\left(k+1\right)\overline{g\left(k-l+1\right)}\]
%   PCONV('rr',f,g) computes the alternative where both f and g are reversed
%   given by
%
%M              L-1
%M     h(l+1) = sum conj(f(-k+1)) * conj(g(k-l+1))
%M              k=0
%     
%F  \[h\left(l+1\right)=\sum_{k=0}^{L-1}f\left(-k+1\right)g\left(l-k+1\right)\]
%   In the above formulas, $l-k$, $k-l$ and $-k$ are computed modulo $L$.
%
%   See also: dft, involute

%   AUTHOR: Peter Soendergaard
%   TESTING: TEST_PCONV
%   REFERENCE: REF_PCONV

% Assert correct input.
error(nargchk(2,3,nargin));

if nargin==2
  ctype='';
  f=p1;
  g=p2;
else
  ctype=p1;
  f=p2;
  g=p3;
end;

if ~all(size(f)==size(g))
  error('f and g must have the same size.');
end;

% The following is a HACK to work around broken support for switch
% statements of empty strings in some Octave versions.
if strcmp(ctype,'')
  ctype=' ';
end;

switch(lower(ctype))
  case {' '}
    h=ifft(fft(f).*fft(g));
  case {'r'}
    h=ifft(fft(f).*conj(fft(g)));
  case {'rr'}
    h=ifft(conj(fft(f)).*conj(fft(g)));
 otherwise
    error('Unknown convolution type: %s',ctype);
end;

% Clean output if input was real-valued
if isreal(f) && isreal(g)
  h=real(h);
end;