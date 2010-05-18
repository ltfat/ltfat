Lr=[9,10,12];

test_failed=0;

disp(' ===============  TEST_HERMBASIS ==========');

for ii=1:length(Lr)
  
  L=Lr(ii);

  H=hermbasis(L);
  
  r1=(H*H')-eye(L);
  res=norm(r1,'fro');
  [test_failed,fail]=ltfatdiditfail(res,test_failed);
  
  s=sprintf('HERMBASIS orth L:%3i %0.5g %s',L,res,fail);
  disp(s);

end;