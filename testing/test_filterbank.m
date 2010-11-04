M=6;
a=3;
N=4;

L=a*N;

g=cell(1,M);
for ii=1:M
  g{ii}=rand(L,1);
end;

gd = filterbankdual(g,a);

f=rand(L,1);

c=filterbank(f,g,a);
r=ifilterbank(c,gd,a);

norm(f-r)

[AF,BF]=filterbankbounds(g,a);

AF
BF

gt = filterbanktight(g,a);

[AF,BF]=filterbankbounds(gt,a);

AF
BF
