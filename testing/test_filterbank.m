M=6;
a=3;
N=4;

L=a*N;

g=cell(1,M);
for ii=1:M
  g{ii}=crand(L,1);
end;

gd = filterbankdual(g,a);

f=rand(L,1);

c_u  = ufilterbank(f,g,a);
c_nu = filterbank(f,g,a);

res=0;
for m=1:M
  res=res+norm(c_nu{m}-c_u(:,m));  
end;

res
r=iufilterbank(c_u,gd,a);

norm(f-r)

[AF,BF]=filterbankbounds(g,a);

AF
BF

gt = filterbanktight(g,a);

[AF,BF]=filterbankbounds(gt,a);

AF
BF

gdreal=filterbankrealdual(g,a);
gtreal=filterbankrealtight(g,a);

rreal=2*real(iufilterbank(c_u,gdreal,a));
[AF,BF]=filterbankrealbounds(gtreal,a);

AF
BF

norm(f-rreal)

ct     = ufilterbank(f,gtreal,a);
rrealt = 2*real(iufilterbank(ct,gtreal,a));

norm(f-rrealt)
