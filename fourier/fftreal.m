function f=fftreal(f,N,dim);
%FFTREAL   FFT for real valued input data
%   Usage: f=fftreal(f);
%          f=fftreal(f,N,dim);
%
%   `fftreal(f)` computes the coefficients corresponding to the positive
%   frequencies of the FFT of the real valued input signal *f*.
%   
%   The function takes exactly the same arguments as fft. See the help on
%   fft for a thorough description.
%
%   See also:  ifftreal, dft

%   AUTHOR    : Peter L. Søndergaard
%   TESTING   : TEST_PUREFREQ
%   REFERENCE : OK
  
complainif_argnonotinrange(nargin,1,3,mfilename);

if nargin<3
  dim=[];  
end;

if nargin<2
  N=[];
end;

if ~isreal(f)
  error('Input signal must be real.');
end;


[f,N,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,N,dim,'FFTREAL');

if ~isempty(N)
   f=postpad(f,N);
end

N2=floor(N/2)+1;

f=comp_fftreal(f);

% Set the new size in the first dimension.
permutedsize(1)=N2;

f=assert_sigreshape_post(f,dim,permutedsize,order);

