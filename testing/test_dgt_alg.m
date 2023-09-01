function test_failed=test_dgt_alg

test_failed=0;

disp(' ===============  TEST_ALG ================');

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
  
  f=tester_crand(L,1);
  g=tester_crand(L,1);
  
  %gd=gabdual(g,a,M);
  
  cc  = ref_dgt(f,g,a,M);
  
  for jj=1:6
    
    cc_comp = feval(['ref_dgt_',num2str(jj)],f,g,a,M);
    
    cdiff=cc-cc_comp;

    res=norm(cdiff(:));      

    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    
    s=sprintf('REF%s L:%3i a:%3i b:%3i c:%3i d:%3i p:%3i q:%3i   %0.5g',...
        num2str(jj),L, a,b,c,d,p,q,res);
    
    disp(s)
    
  end;
end;




