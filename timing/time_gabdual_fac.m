a=300;
M=400;

L=a*M;

N=L/a;
b=L/M;

c=gcd(a,M);
p=a/c;
q=M/c;
d=N/q;

gf=crand(p*q,c*d);

tic;
gdf1=comp_gabdual_fac(gf,L,a,M);
toc

tic
gdf2=ref_gabdual_fac_time(gf,L,a,M);
toc


