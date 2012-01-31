function c=ref_wmdct2(f,g1,g2,M1,M2);

L1=size(f,1);
L2=size(f,2);

c=wmdct(f,g1,M1);

c=reshape(c,L1,L2);

c=c.';

c=wmdct(c,g2,M2);

c=reshape(c,L2,L1);

c=c.';

c=reshape(c,M1,L1/M1,M2,L2/M2);

