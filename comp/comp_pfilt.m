function h=comp_pfilt(f,g,a)
%COMP_PFILT  Compute PFILT

[L W]=size(f);
  
g=fir2long(g,L);
N=L/a;
L2=floor(L/2)+1;
N2=floor(N/2)+1;
a2=floor(a/2)+1;

% Force FFT along dimension 1, since we have permuted the dimensions
% manually
if a==1
  if isreal(f) && isreal(g)
    h=ifftreal(fftreal(f,L,1).*repmat(fftreal(g,L,1),1,W),L,1);
  else
    h=ifft(fft(f,L,1).*repmat(fft(g,L,1),1,W),L,1);
  end;
else
  % Use Poisson

  if isreal(f) && isreal(g)
    % Algorithm can be optimized, current method computes the full spectrum by
    % FFT, and then only extracts the desired part after summation.
    tmp=fft(f,L,1).*repmat(fft(g,L,1),1,W);
    tmp=squeeze(sum(reshape(tmp,N,a,W),2));
           
    h=ifftreal(tmp(1:N2,:),N)/a;
    
  else
    h=squeeze(ifft(sum(reshape(fft(f,L,1).*repmat(fft(g,L,1),1,W),N,a,W),2)))/a;
  end;

  
end;

