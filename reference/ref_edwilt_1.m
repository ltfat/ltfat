function c=ref_edwilt_1(f,g,M)
%REF_EDWILT_1  Reference Even DWILT by DWILTII
%
%  The coefficients cannot be arranged in a rectangular layout.
L=size(g,1)/2;
W=size(f,2);

a=M;
N=L/a;

f2=[f;flipud(conj(f))];

c2=ref_dwiltii(f2,g,M,M);

reshape(c2,2*M,N)

c=zeros(M*N,1);

c(1:

