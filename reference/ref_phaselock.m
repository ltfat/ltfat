function c = ref_phaselock(c,a)
%REF_PHASELOCK  Phaselock Gabor coefficients
%   Usage:  c=phaselock(c,a);
%
%   `phaselock(c,a)` phaselocks the Gabor coefficients *c*. The coefficients
%   must have been obtained from a |dgt| with parameter *a*.

%   AUTHOR:

M=size(c,1);
N=size(c,2);
L=N*a;
b=L/M;

TimeInd = (0:(N-1))*a;
FreqInd = (0:(M-1))*b;

phase = FreqInd'*TimeInd;
phase = exp(2*1i*pi*phase/L);

c=bsxfun(@times,c,phase);

