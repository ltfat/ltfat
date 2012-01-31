function c=ref_dgt_pconv(f,g,a,M);
%REF_DGT_PCONV  Compute a DGT using PCONV

  
L=size(f,1);

N=L/a;
b=L/M;

c=zeros(M,N);

for m=0:M-1
  work = pconv(f,involute(g).*expwave(L,m*b));
  c(m+1,:) = reshape(work(1:a:L),1,N);
end;

c=phaseunlock(c,a);

