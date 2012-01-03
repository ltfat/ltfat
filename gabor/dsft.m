function F=dsft(F);
%DSFT  Discrete Symplectic Fourier Transform
%   Usage:  C=dsft(F);
%
%   DSFT(F) computes the discrete symplectic Fourier transform of F.
%   F must be a matrix or a 3D array. If F is a 3D array, the 
%   transformation is applied along the first two dimensions.
%
%   Let F be a _K x _L matrix. Then the DSFT of F is given by
%
%M                              L-1 K-1
%M   C(m+1,n+1) = 1/sqrt(K*L) * sum sum F(k+1,l+1)*exp(2*pi*i(k*n/K-l*m/L))
%M                              l=0 k=0
%
%F \[C\left(m+1,n+1\right)=\frac{1}{\sqrt{KL}}\sum_{l=0}^{L-1}\sum_{k=0}^{K-1}F
%F \left(k+1,l+1\right)e^{2\pi i\left(kn/K-lm/L\right)}\]
%   for $m=0,...,L-1$ and $n=0,...,K-1$.
%
%   The DSFT is its own inverse.
%
%   References: fekolu06 

error(nargchk(1,1,nargin));

D=ndims(F);

if (D<2) || (D>3)
  error('Input must be two/three dimensional.');
end;

W=size(F,3);

if W==1
  F=dft(idft(F).');
else
  % Apply to set of planes.
  
  R1=size(F,1);
  R2=size(F,2);
  Fo=zeros(R2,R1,W);
  for w=1:W
    Fo(:,:,w)=dft(idft(F(:,:,w).'));
  end;
  F=Fo;
end;
%OLDFORMAT
