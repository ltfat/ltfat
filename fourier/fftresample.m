function f=fftresample(f,L,dim)
%FFTRESAMPLE   Resample signal using Fourier interpolation
%   Usage:  h=fftresample(f,L);
%           h=fftresample(f,L,dim);
%
%   FFTRESAMPLE(f,L) returns a Fourier interpolation of the signal f
%   to length L. If the function is applied to a matrix, it will apply
%   to each column.  
%
%   FFTRESAMPLE(f,L,dim) does the same along dimension dim.
%
%   If the input signal is NOT a periodic signal (or close to), the
%   DCTRESAMPLE method will give much better results at the
%   endpoints.
%
%   See also:  dctresample, middlepad

%   AUTHOR: Peter Soendergaard
  
% ------- Checking of input --------------------
error(nargchk(2,3,nargin));

if nargin<3
  dim=[];
end;

[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,L,dim,'FFTRESAMPLE');

wasreal=isreal(f);

% The 'dim=1' below have been added to avoid fft and middlepad being
% smart about choosing the dimension.
f=ifft(middlepad(fft(f,[],1),L,1))/Ls*L;

f=assert_sigreshape_post(f,dim,permutedsize,order);

if wasreal
  f=real(f);
end;
%OLDFORMAT
