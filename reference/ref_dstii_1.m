function c=ref_dstii_1(f)
%REF_DSTII  Reference Discrete Sine Transform type II
%   Usage:  c=ref_dstii(f);
%
%   The transform is computed by an FFT of 4 times the length of f.
%

L=size(f,1);
W=size(f,2);

if ~isreal(f)
  c=ref_dstii_1(real(f)) + i*ref_dstii_1(imag(f));
else
  lf=zeros(L*4,W);
  
  lf(2:2:2*L,:)=f;
  lf(2*L+2:2:end,:)=-flipud(f);

  fflong=imag(fft(lf));

  c=-fflong(2:L+1,:)/sqrt(2*L);
  
  % Scale last coefficients to obtain orthonormal transform.
  c(end,:)=1/sqrt(2)*c(end,:);
end

