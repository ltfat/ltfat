function c=ref_dctiii_1(f)
%REF_DCTII  Reference Discrete Consine Transform type III
%   Usage:  c=ref_dctiii(f);
%
%   The transform is computed by an FFT of 4 times the length of f.
%   See the Wikipedia article on "Discrete Cosine Transform"
%
%   This is the inverse of REF_DCTII

L=size(f,1);
W=size(f,2);

if ~isreal(f)
  c=ref_dctiii_1(real(f))+i*ref_dctiii_1(imag(f));
else
  % Scale coefficients to obtain orthonormal transform.
  f(1,:)=sqrt(2)*f(1,:);
  f=f*sqrt(2*L);
  
  % Make 4x long vector
  lf=[f;...
      zeros(1,W);...
      -flipud(f);...
      -f(2:end,:);...
      zeros(1,W);...
      flipud(f(2:end,:))];
  

  fflong=real(ifft(lf));

  c=fflong(2:2:2*L,:);
end

