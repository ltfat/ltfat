function test_failed=test_involute

Lr=[9,10];

test_failed=0;

disp(' ===============  TEST_INVOLUTE ===========');

for ii=1:length(Lr)
  
  L=Lr(ii);
  f=tester_crand(L,1);
  
  r1=conj(dft(f));
  r2=dft(involute(f));
  
  res=norm(r1-r2);
  [test_failed,fail]=ltfatdiditfail(res,test_failed);          
  s=sprintf('INVOLUTE  L:%3i %0.5g %s',L,res,fail);
  disp(s);

end;

