function F = comp_ifilterbank_fft(c,G,a)
%COMP_IFILTERBANK_FFT  Compute filtering in FD
%
%

W = size(c{1},2);
M = numel(G);
L = numel(G{1});
F = zeros(L,W,assert_classname(c{1},G{1}));

for m=1:M
   for w=1:W
      % This repmat cannot be replaced by bsxfun
      F(:,w)=F(:,w)+repmat(fft(c{m}(:,w)),a(m),1).*conj(G{m});
   end;
end
