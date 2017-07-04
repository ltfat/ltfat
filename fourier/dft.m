function f=dft(f,N,dim);
%DFT   Normalized Discrete Fourier Transform
%   Usage: f=dft(f);
%          f=dft(f,N,dim);
%
%   `dft` computes a normalized or unitary discrete Fourier transform. The 
%   unitary discrete Fourier transform is computed by
%   
%   ..                     L-1
%     c(k+1) = 1/sqrt(L) * sum f(l+1)*exp(-2*pi*i*k*l/L)
%                          l=0
%
%   .. math:: c\left(k+1\right)=\frac{1}{\sqrt{L}}\sum_{l=0}^{L-1}f\left(l+1\right)e^{-2\pi ikl/L}
%
%   for $k=0,\ldots,L-1$.
%
%   The output of `dft` is a scaled version of the output from `fft`. The
%   function takes exactly the same arguments as `fft`. See the help on `fft`
%   for a thorough description.
%
%   See also:  idft

%   AUTHOR: Peter L. Søndergaard, Jordy van Velthoven
%   TESTING: TEST_DFT
%   REFERENCE: REF_DFT

complainif_argnonotinrange(nargin,1,3,mfilename);

if nargin<3
  dim=[];
end;

if nargin<2
  N=[];
end;

[f,N,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,N,dim,'DFT');

% Force FFT along dimension 1, since we have permuted the dimensions
% manually
f=fft(f,N,1)/sqrt(N);

f=assert_sigreshape_post(f,dim,permutedsize,order);

