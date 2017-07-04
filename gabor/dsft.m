function F=dsft(F);
%DSFT  Discrete Symplectic Fourier Transform
%   Usage:  C=dsft(F);
%
%   `dsft(F)` computes the discrete symplectic Fourier transform of *F*.
%   *F* must be a matrix or a 3D array. If *F* is a 3D array, the 
%   transformation is applied along the first two dimensions.
%
%   Let *F* be a $K \times L$ matrix. Then the DSFT of *F* is given by
%
%   ..                           L-1 K-1
%     C(m+1,n+1) = 1/sqrt(K*L) * sum sum F(k+1,l+1)*exp(2*pi*i(k*n/K-l*m/L))
%                                l=0 k=0
%
%   .. math:: C\left(m+1,n+1\right)=\frac{1}{\sqrt{KL}}\sum_{l=0}^{L-1}\sum_{k=0}^{K-1}F
%            \left(k+1,l+1\right)e^{2\pi i\left(kn/K-lm/L\right)}
%
%   for $m=0,\ldots,L-1$ and $n=0,\ldots,K-1$.
%
%   The `dsft` is its own inverse.
%
%   References: feichtinger2008metaplectic

%   AUTHOR : Peter L. Søndergaard, Jordy van Velthoven (TESTING).
%   TESTING: TEST_DSFT 
%   REFERENCE: REF_DSFT

complainif_argnonotinrange(nargin,1,1,mfilename);

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
  Fo=zeros(R2,R1,W,assert_classname(F));
  for w=1:W
    Fo(:,:,w)=dft(idft(F(:,:,w).'));
  end;
  F=Fo;
end;

