a=20;
M=30;

L=a*M;
N=L/a;

c=randn(M,N);

g=gabtight(a,M,L);

[V,D]=gabmuleigs(10,c,g,a);