function h=lconv(f,g,varargin)
%LCONV  Linear convolution
%   Usage:  h=lconv(f,g);
%
%   `lconv(f,g)` computes the linear convolution of *f* and *g*. The linear 
%   convolution is given by
%
%   ..          Lh-1
%      h(l+1) = sum f(k+1) * g(l-k+1)
%               k=0
%
%   .. math:: h\left(l+1\right)=\sum_{k=0}^{L_{h}-1}f\left(k+1\right)g\left(l-k+1\right)
%
%   with $L_{h} = L_{f} + L_{g} - 1$ where $L_{f}$ and $L_{g}$ are the lengths of *f* and *g*, 
%   respectively.
%
%   `lconv(f,g,'r')` computes the linear convolution of *f* and *g* where *g* is reversed.
%   This type of convolution is also known as linear cross-correlation and is given by
%
%   ..          Lh-1
%      h(l+1) = sum f(k+1) * conj(g(k-l+1))
%               k=0
%
%   .. math:: h\left(l+1\right)=\sum_{k=0}^{L_{h}-1}f\left(k+1\right)\overline{g\left(k-l+1\right)}
%
%   `lconv(f,g,'rr')` computes the alternative where both *f* and *g* are
%   reversed given by
%
%   ..          Lh-1
%      h(l+1) = sum conj(f(-k+1)) * conj(g(k-l+1))
%               k=0
%     
%   .. math:: h\left(l+1\right)=\sum_{k=0}^{L_{h}-1}f\left(-k+1\right)g\left(l-k+1\right)
%
%   In the above formulas, $l-k$, $k-l$ and $-k$ are computed modulo $L_{h}$.
%
%   The input arrays *f* and *g* can be 1D vectors or one of them can be
%   a multidimensional array. In either case, the convolution is performed
%   along columns with row vectors transformed to columns.
%
%   See also: pconv

%   AUTHOR: Jordy van Velthoven
%   TESTING: TEST_LCONV	
%   REFERENCE: REF_LCONV
  
complainif_notenoughargs(nargin, 2, 'LCONV');

definput.keyvals.L=[];
definput.keyvals.dim=[];
definput.flags.type={'default', 'r', 'rr'};

[flags,~,L,dim]=ltfatarghelper({'L','dim'},definput,varargin);

[f,~,Lf,Wf,dimoutf,permutedsize_f,order_f]=assert_sigreshape_pre(f,L,dim,'LCONV');
[g,~,Lg,Wg,dimoutg,permutedsize_g,order_g]=assert_sigreshape_pre(g,L,dim,'LCONV');

if (Wf>1) && (Wg>1)
  error('%s: Only one of the inputs can be multi-dimensional.',upper(mfilename));
end;

W=max(Wf,Wg);
if Wf<W
  f=repmat(f,1,W);
end;

if Wg<W
  g=repmat(g,1,W);
end;

Lh = Lf+Lg-1;

f = postpad(f,Lh);
g = postpad(g,Lh);

if isreal(f) && isreal(g)
  fftfunc = @(x) fftreal(x);
  ifftfunc = @(x) ifftreal(x, Lh);
else
  fftfunc = @(x) fft(x);
  ifftfunc = @(x) ifft(x, Lh);
end;

if flags.do_default
  h=ifftfunc(fftfunc(f).*fftfunc(g));
end;

if flags.do_r
  h=ifftfunc(fftfunc(f).*(conj(fftfunc(g))));
end;

if flags.do_rr
  h=ifftfunc((conj(fftfunc(f))).*(conj(fftfunc(g))));
end;


