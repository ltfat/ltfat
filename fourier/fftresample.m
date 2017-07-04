function f=fftresample(f,L,dim)
%FFTRESAMPLE   Resample signal using Fourier interpolation
%   Usage:  h=fftresample(f,L);
%           h=fftresample(f,L,dim);
%
%   `fftresample(f,L)` returns a Fourier interpolation of the signal *f*
%   to length *L*. If the function is applied to a matrix, it will apply
%   to each column.  
%
%   `fftresample(f,L,dim)` does the same along dimension *dim*.
%
%   If the input signal is **not** a periodic signal (or close to), the
%   |dctresample| method gives much better results at the endpoints.
%
%   See also:  dctresample, middlepad

%   AUTHOR: Peter L. Søndergaard
  
% ------- Checking of input --------------------
complainif_argnonotinrange(nargin,2,3,mfilename);

if nargin<3
  dim=[];
end;

[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,L,dim,'FFTRESAMPLE');

wasreal=isreal(f);

% The 'dim=1' below have been added to avoid fft and middlepad being
% smart about choosing the dimension.
% In addition, postpad is explicitly told to pad with zeros.
if isreal(f)
  L2=floor(L/2)+1;
  f=ifftreal(postpad(fftreal(f,[],1),L2,0,1),L,1)/Ls*L;
else
  f=ifft(middlepad(fft(f,[],1),L,1))/Ls*L;
end;

f=assert_sigreshape_post(f,dim,permutedsize,order);

