function c=ref_dctii_1(f)
%REF_DCTII  Reference Discrete Consine Transform type II
%   Usage:  c=ref_dctii(f);
%
%   The transform is computed by an FFT of 4 times the length of f.
%   See the Wikipedia article on "Discrete Cosine Transform"
%

L=size(f,1);
W=size(f,2);

if ~isreal(f)
  c=ref_dctii_1(real(f))+i*ref_dctii_1(imag(f));
else
  lf=zeros(L*4,W);
  
  lf(2:2:2*L,:)=f;
  lf(2*L+2:2:end,:)=flipud(f);

  fflong=real(fft(lf));

  c=fflong(1:L,:)/sqrt(2*L);
  
  % Scale first coefficients to obtain orthonormal transform.
  c(1,:)=1/sqrt(2)*c(1,:);
end;

