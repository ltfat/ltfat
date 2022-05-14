function c=comp_filterbank_fft(F,G,a)
%COMP_FILTERBANK_FFT  Compute filtering in FD
%
%   does the same as comp_filterbank_fftbl, but for
%   filters that are not bandlimited in the frequency domain

M = numel(G);
[L,W] = size(F);
c = cell(M,1);
N = L./a;

for m=1:M
    c{m}=zeros(N(m),W,assert_classname(F,G{m}));
    for w=1:W
        c{m}(:,w)=ifft(sum(reshape(F(:,w).*G{m},N(m),a(m)),2))/a(m);
    end;
end
