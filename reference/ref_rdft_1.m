function c=ref_rdft_1(f)
%REF_RDFT_1  Reference RDFT by FFT
%   Usage:  c=ref_rdft_1(f);
%
%   Compute RDFT by doing a DFT and returning half the coefficients.
%
%   The transform is orthonormal



L=size(f,1);
Lhalf=ceil(L/2);
Lend=Lhalf*2-1;

if ~isreal(f)
  c=ref_rdft_1(real(f))+i*ref_rdft_1(imag(f));
else
  cc=fft(f);
  
  c=zeros(size(f));
  
  % Copy the first coefficient, it is real
  c(1,:)=1/sqrt(2)*real(cc(1,:));
  
  % Copy the cosine-part of the coefficients.
  c(2:2:Lend,:)=real(cc(2:Lhalf,:));
  
  % Copy the sine-part of the coefficients.
  c(3:2:Lend,:)=-imag(cc(2:Lhalf,:));
  
  % If f has an even length, we must also copy the Niquest-wave
  % (it is real)
  if mod(L,2)==0
    c(end,:)=1/sqrt(2)*real(cc(L/2+1,:));
  end;
  
  % Make it an ortonomal transform
  c=c/sqrt(L/2);
end;

