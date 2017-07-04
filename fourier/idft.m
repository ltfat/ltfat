function f=idft(c,N,dim)
%IDFT  Inverse normalized Discrete Fourier Transform
%   Usage: f=idft(c);
%          f=idft(c,N,dim);
%
%   `idft` computes a normalized or unitary inverse discrete Fourier transform. 
%   The unitary discrete Fourier transform is computed by
%   
%   ..                     L-1
%     f(l+1) = 1/sqrt(L) * sum c(k+1)*exp(2*pi*i*k*l/L)
%                          k=0
%
%   .. math:: f\left(l+1\right)=\frac{1}{\sqrt{L}}\sum_{k=0}^{L-1}c\left(k+1\right)e^{2\pi ikl/L}
%
%   for $l=0,\ldots,L-1$.
%
%   The output of `idft` is a scaled version of the output from `ifft`. The
%   function takes exactly the same arguments as `ifft`. See the help on `ifft`
%   for a thorough description.
%
%   See also:  dft

%   AUTHOR: Peter L. Søndergaard, Jordy van Velthoven
%   TESTING: TEST_IDFT
%   REFERENCE: TEST_DFT

complainif_argnonotinrange(nargin,1,3,mfilename);

if nargin<3
  dim=[];  
end;

if nargin<2
  N=[];
end;

[c,N,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(c,N,dim,'IDFT');

% Force IFFT along dimension 1, since we have permuted the dimensions
% manually
f=ifft(c,N,1)*sqrt(N);

f=assert_sigreshape_post(f,dim,permutedsize,order);

