function f=ifftreal(c,N,dim);
%IFFTREAL  Inverse FFT for real valued signals
%   Usage: f=ifftreal(c,N);
%          f=ifftreal(c,N,dim);
%
%   `ifftreal(c,N)` computes an inverse FFT of the positive frequency
%   Fourier coefficients *c*. The length *N* must always be specified,
%   because the correct transform length cannot be determined from the
%   size of *c*.
%
%   `ifftreal(c,N,dim)` does the same along dimension *dim*.
%
%   See also:  fftreal

%   AUTHOR : Peter L. Søndergaard

complainif_argnonotinrange(nargin,2,3,mfilename);

if nargin==2
  dim=[];  
end;

N2=floor(N/2)+1;

[c,~,~,~,dim,permutedsize,order]=assert_sigreshape_pre(c,N2,dim,'IFFTREAL');

% Clean for safety
c(1,:)=real(c(1,:));

f=comp_ifftreal(c,N);

% Restore the full size in the first dimension.
permutedsize(1)=N;

f=assert_sigreshape_post(f,dim,permutedsize,order);
