function f=dft(f,N,dim);
%DFT   Normalized Discrete Fourier Transform
%   Usage: f=dft(f);
%          f=dft(f,N,dim);
%
%   DFT computes a normalized discrete Fourier transform. This is nothing
%   but a scaled version of the output from FFT. The function takes exactly
%   the same arguments as FFT. See the help on FFT for a thorough
%   description.
%
%   See also:  idft

%   AUTHOR: Peter Soendergaard
%   TESTING: OK
%   REFERENCE: OK

error(nargchk(1,3,nargin));

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
