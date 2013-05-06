function c=ref_dft_1(f)
%REF_DFT_1  Reference DFT by doubling 
%   Usage:  c=ref_dft_1(f);
%
%   REF_DFT_1(f) computes a DFT of f by upsampling f and inserting zeros
%   at the odd positions.
%
%   This is not an efficient method, it is just meant to illustrate a 
%   symmetry of the DFT.

L=size(f,1);
W=size(f,2);

flong=zeros(2*L,W,assert_classname(f));
flong(1:2:end-1)=f;

fflong=fft(flong)/sqrt(L);

c=fflong(1:L,:);


