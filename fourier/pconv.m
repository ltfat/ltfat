function h=pconv(f,g,varargin)
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
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~all(size(f)==size(g))
  error('f and g must have the same size.');
end;

definput.flags.type={'default','r','rr'};

[flags,kv]=ltfatarghelper({},definput,varargin);

if flags.do_default
    h=ifft(fft(f).*fft(g));
end;

if flags.do_r
  h=ifft(fft(f).*conj(fft(g)));
end;

if flags.do_rr
  h=ifft(conj(fft(f)).*conj(fft(g)));
end;

% Clean output if input was real-valued
if isreal(f) && isreal(g)
  h=real(h);
end;