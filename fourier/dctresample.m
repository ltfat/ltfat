function f=dctresample(f,L,dim)
%DCTRESAMPLE   Resample signal using Fourier interpolation
%   Usage:  h=dctresample(f,L);
%           h=dctresample(f,L,dim);
%
%   `dctresample(f,L)` returns a discrete cosine interpolation of the signal *f*
%   to length *L*. If the function is applied to a matrix, it will apply
%   to each column.
%
%   `dctresample(f,L,dim)` does the same along dimension *dim*.
%
%   If the input signal is not a periodic signal (or close to), this method
%   will give much better results than |fftresample| at the endpoints, as
%   this method assumes than the signal is even a the endpoints.
%
%   The algorithm uses a DCT type iii.
%
%   See also:  fftresample, middlepad, dctiii

%   AUTHOR: Peter L. Søndergaard
  
% ------- Checking of input --------------------
complainif_argnonotinrange(nargin,2,3,mfilename);

if nargin<3
  dim=[];
end;

[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,L,dim,'DCTRESAMPLE');

wasreal=isreal(f);

% The 'dim=1' below have been added to avoid dct and middlepad being
% smart about choosing the dimension.
f=dctiii(postpad(dctii(f,[],1),L))*sqrt(L/Ls);

f=assert_sigreshape_post(f,dim,permutedsize,order);

if wasreal
  f=real(f);
end;

