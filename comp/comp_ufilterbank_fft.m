function c=comp_ufilterbank_fft(f,g,a);  
%COMP_UFILTERBANK_FFT   Classic filtering by FFT
%   Usage:  c=comp_ufilterbank_fft(f,g,a);
%

L=size(f,1);
W=size(f,2);
M=size(g,2);

N=L/a;

c=zeros(N,M,W,assert_classname(f,g));

% This routine does not yet use FFTREAL, because it must be able to
% handle downsampling, which is much easier to express in the FFT case.
G=fft(fir2long(g,L));

for w=1:W
  F=fft(f(:,w));
  for m=1:M
    c(:,m,w)=ifft(sum(reshape(F.*G(:,m),N,a),2))/a;
  end;
end;

if isreal(f) && isreal(g)
  c=real(c);
end;
  

