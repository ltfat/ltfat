function c=ref_dstiii_1(f)
%REF_DSTII  Reference Discrete Sine Transform type III
%   Usage:  c=ref_dstiii(f);
%
%   This is the inverse of REF_DSTII

L=size(f,1);
W=size(f,2);

if ~isreal(f)
  c=ref_dstiii_1(real(f))+i*ref_dstiii_1(imag(f));
else

  % Scale coefficients to obtain orthonormal transform.
  f(end,:)=sqrt(2)*f(end,:);
  f=-f*sqrt(2*L);

  % Make 4x long vector
  lf=[zeros(1,W);...
      i*f;...
      i*flipud(f(1:end-1,:));...
      zeros(1,W);...
      -i*f;...
      -i*flipud(f(1:end-1,:));...
      ];
  
  fflong=real(ifft(lf));

  c=fflong(2:2:2*L,:);

end

