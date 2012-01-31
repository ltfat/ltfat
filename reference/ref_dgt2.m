function c=ref_dgt2(f,g1,g2,a1,a2,M1,M2);
%REF_DGT2  Reference DGT2
%
%  Compute a DGT2 using a DGT along each dimension.

L1=size(f,1);
L2=size(f,2);

N1=L1/a1;
N2=L2/a2;

c=dgt(f,g1,a1,M1);

c=reshape(c,M1*N1,L2);

c=c.';

c=dgt(c,g2,a2,M2);

c=reshape(c,M2*N2,M1*N1);

c=c.';

c=reshape(c,M1,N1,M2,N2);

