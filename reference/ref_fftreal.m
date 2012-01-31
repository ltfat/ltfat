function c=ref_fftreal(f)
%REF_FFTREAL  Reference FFTREAL
%
%  FFTREAL is computed by doing an FFT, and keeping only the positive
%  frequencies.
  

L=size(f,1);
c=fft(f);

c=c(1:floor(L/2)+1,:);


