function f=comp_fftreal(f)
%COMP_FFTREAL  Compute an FFTREAL
  
N=size(f,1);
N2=floor(N/2)+1;

% Force FFT along dimension 1, since we have permuted the dimensions
% manually
f=fft(f,N,1);
f=f(1:N2,:);


