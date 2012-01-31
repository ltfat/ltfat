function c=ref_dfft(f,order)
%REF_DFFT  Reference Discrete Fractional Fourier Transform
%   Usage:  c=ref_dfft(f,order);
%

L=size(f,1);

% Create matrix representation of the DFT
F=idft(eye(L));

c=(F^order)*f;


