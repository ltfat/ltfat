function c=ref_dftiv_1(f)
%REF_DFT_1  Reference DFTIV by quadrupling
%   Usage:  c=ref_dftii_1(f);
%
%   REF_DFTIV_1(f) computes a DFTIV of f by upsampling f and inserting zeros
%   at the even positions, and then doubling this signal

L=size(f,1);
W=size(f,2);

flong=zeros(4*L,W);
flong(2:2:2*L,:)=f;

fflong=fft(flong)/sqrt(L);

c=fflong(2:2:2*L,:);


