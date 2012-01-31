
L=24;
Lg=4;

Lb=8;

f=randn(L,1);
g=randn(Lg,1);

ref1=pconv(f,postpad(g,L));
ola1=ref_pconv_ola_postpad(f,g,Lb);

norm(ref1-ola1)

ref2=pconv(f,fir2long(g,L));
ola2=ref_pconv_ola_fir2long(f,g,Lb);

norm(ref2-ola2)

[ref2,ola2,ref2-ola2]

