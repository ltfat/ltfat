function c=ref_dwilt2(f,g1,g2,M1,M2);

L1=size(f,1);
L2=size(f,2);

c=dwilt(f,g1,M1);

c=reshape(c,L1,L2);

c=c.';

c=dwilt(c,g2,M2);

c=reshape(c,L2,L1);

c=c.';

c=reshape(c,M1*2,L1/M1/2,M2*2,L2/M2/2);

