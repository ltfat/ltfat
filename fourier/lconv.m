function h=lconv(f,g,varargin)
%LCONV  Linear convolution
%   Usage:  h=lconv(f,g);
%
%   `lconv(f,g)` computes the linear convolution of *f* and *g*. The linear 
%   convolution is given by
%
%   ..          L-1
%      h(l+1) = sum f(k+1) * g(l-k)+1)
%               k=0
%
%   .. math:: h\left(l+1\right)=\sum_{k=0}^{L-1}f\left(k+1\right)g\left(l-k+1\right)
%
%   `lconv(f,g,'r')` computes the linear convolution of *f* and *g* where *g* is reversed.
%   This type of convolution is also known as linear cross-correlation and is given by
%
%   ..          L-1
%      h(l+1) = sum f(k+1) * conj(g(k-l+1))
%               k=0
%
%   .. math:: h\left(l+1\right)=\sum_{k=0}^{L-1}f\left(k+1\right)\overline{g\left(k-l+1\right)}
%
%   `lconv(f,g,'rr')` computes the alternative where both *f* and *g* are
%   reversed given by
%
%   ..          L-1
%      h(l+1) = sum conj(f(-k+1)) * conj(g(k-l+1))
%               k=0
%     
%   .. math:: h\left(l+1\right)=\sum_{k=0}^{L-1}f\left(-k+1\right)g\left(l-k+1\right)
%
%   Since all these types of convolution are linear, the output of `lconv` will have a length L ::
%
%     L = length(f)+length(g)-1;
%
%   See also: pconv

%   AUTHOR: Jordy van Velthoven
%	TESTING: TEST_LCONV	
%	REFERENCE: REF_LCONV
  
complainif_notenoughargs(nargin, 2, 'LCONV');

definput.keyvals.L=[];
definput.keyvals.dim=[];
definput.flags.type={'default', 'r', 'rr'};

[flags,kv,L,dim]=ltfatarghelper({'L','dim'},definput,varargin);

[f,L1,Lf,Wf,dimout,permutedsize_f,order_f]=assert_sigreshape_pre(f,L,dim,'LCONV');
[g,L2,Lg,Wg,dimout,permutedsize_g,order_g]=assert_sigreshape_pre(g,L,dim,'LCONV');

Lh = Lf+Lg-1;

f = [f; zeros(Lh - Lf, 1)];
g = [g; zeros(Lh - Lg, 1)];

if isreal(f) && isreal(g)
  fftfunc = @(x) fftreal(x);
  ifftfunc = @(x) ifftreal(x, Lh);
else
  fftfunc = @(x) fft(x);
  ifftfunc = @(x) ifft(x, Lh);
end;

if flags.do_default
  h=ifftfunc(fftfunc(f).*fftfunc(g), Lh);
end;

if flags.do_r
  h=ifftfunc(fftfunc(f).*(conj(fftfunc(g))), Lh);
end;

if flags.do_rr
  h=ifftfunc((conj(fftfunc(f))).*(conj(fftfunc(g))), Lh);
end;
