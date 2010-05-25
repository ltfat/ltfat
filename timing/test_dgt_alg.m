Lr=[24, 24, 24, 144,108,144,24,135,35,77,20];
ar=[ 4,  3,  6,   9,  9, 12, 6,  9, 5, 7, 1];
Mr=[ 6,  4,  4,  16, 12, 24, 8,  9, 7,11,20];

for ii=1:length(Lr);

  L=Lr(ii);
  
  M=Mr(ii);
  a=ar(ii);
  
  b=L/M;
  N=L/a;
  c=gcd(a,M);
  d=gcd(b,N);
  p=a/c;
  q=M/c;
  
  f=crand(L,1);
  g=crand(L,1);
  
  %gd=gabdual(g,a,M);
  
  cc  = ref_dgt(f,g,a,M);
  cc1 = ref_dgt_1(f,g,a,M);
  cc2 = ref_dgt_2(f,g,a,M);
  cc3 = ref_dgt_3(f,g,a,M);
  cc4 = ref_dgt_4(f,g,a,M);
  cc5 = ref_dgt_5(f,g,a,M);
  cc6 = dgt(f,g,a,M);  

  cdiff1=cc-cc1;
  cdiff2=cc-cc2;
  cdiff3=cc-cc3;
  cdiff4=cc-cc4;
  cdiff5=cc-cc5;
  cdiff6=cc-cc6;

  res1=norm(cdiff1(:));      
  res2=norm(cdiff2(:));      
  res3=norm(cdiff3(:));      
  res4=norm(cdiff4(:));      
  res5=norm(cdiff5(:));      
  res6=norm(cdiff6(:));      

  s1=sprintf('REF1 L:%3i a:%3i b:%3i c:%3i d:%3i p:%3i q:%3i   %0.5g',L, ...
             a,b,c,d,p,q,res1);
  s2=sprintf('REF2 L:%3i a:%3i b:%3i c:%3i d:%3i p:%3i q:%3i   %0.5g',L, ...
             a,b,c,d,p,q,res2);
  s3=sprintf('REF3 L:%3i a:%3i b:%3i c:%3i d:%3i p:%3i q:%3i   %0.5g',L, ...
             a,b,c,d,p,q,res3);
  s4=sprintf('REF4 L:%3i a:%3i b:%3i c:%3i d:%3i p:%3i q:%3i   %0.5g',L, ...
             a,b,c,d,p,q,res4);
  s5=sprintf('DGT5 L:%3i a:%3i b:%3i c:%3i d:%3i p:%3i q:%3i   %0.5g',L, ...
             a,b,c,d,p,q,res5);
  s6=sprintf('DGT  L:%3i a:%3i b:%3i c:%3i d:%3i p:%3i q:%3i   %0.5g',L, ...
             a,b,c,d,p,q,res6);

  disp(s1)
  disp(s2)
  disp(s3)
  disp(s4)
  disp(s5)
  disp(s6)
      
end;


