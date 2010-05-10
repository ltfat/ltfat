function f=ifftreal(c,N,dim);
%IFFTREAL  Inverse FFT for real valued signals.
%   Usage: f=ifftreal(c,N);
%          f=ifftreal(c,N,dim);
%
%   IFFTREAL(c,N) computes an inverse FFT of the positive frequency
%   Fourier coefficients c. The length N must always be specified,
%   because the correct transform length cannot be determined from the
%   size of c.
%
%   IFFTREAL(c,N,dim) does the same along dimension dim.
%
%   See also:  fftreal

%   AUTHOR : Peter Soendergaard

error(nargchk(2,3,nargin));

if nargin==2
  dim=[];  
end;

N2=floor(N/2)+1;

[c,N2,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(c,N2,dim,'IFFTREAL');

f=comp_ifftreal(c,N);

% Restore the full size in the first dimension.
permutedsize(1)=N;

f=assert_sigreshape_post(f,dim,permutedsize,order);
