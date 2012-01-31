function c=ref_rdftiii_1(f)
%REF_RDFTIII_1  Reference RDFT by FFT
%   Usage:  c=ref_rdftiii_1(f);
%
%   Compute RDFTII by doing a DFTIII and returning half the coefficients.
%   Only works for real functions.
%
%   The transform is orthonormal


L=size(f,1);
Lhalf=floor(L/2);
Lend=Lhalf*2;

cc=ref_dftiii(f);

c=zeros(size(f));

% Copy the cosine-part of the coefficients.
c(1:2:Lend,:)=sqrt(2)*real(cc(1:Lhalf,:));

% Copy the sine-part of the coefficients.
c(2:2:Lend,:)=-sqrt(2)*imag(cc(1:Lhalf,:));

% If f has an odd length, we must also copy the Niquest-wave
% (it is real)
if mod(L,2)==1
  c(end,:)=real(cc((L+1)/2,:));
end;


