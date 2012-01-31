function c=ref_dctiii_2(f)
%DCTII  Reference Discrete Consine Transform type III
%   Usage:  c=ref_dctiii_2(f);
%
%   The transform is real (only works for real input data) and
%   it is orthonormal.
%
%   The transform is computed as the exact inverse of DCTII, i.e. all
%   steps in the DCTII are reversed in order of computation.
%
%   NOT WORKING

L=size(f,1);
W=size(f,2);

R=1/sqrt(2)*[diag(exp((0:L-1)*pi*i/(2*L)));...
	     zeros(1,L); ...
	     [zeros(L-1,1),flipud(diag(exp(-(1:L-1)*pi*i/(2*L))))]];

R

R(1,1)=1;

c=real(R'*fft([f;flipud(f)])/sqrt(L)/2);


