function test_failed=test_pconv
Lr=[9,10];

ctypes={'default','r','rr'};

test_failed=0;

disp(' ===============  TEST_PCONV ==============');

for jj=1:length(Lr)
  L=Lr(jj);
  for ii=1:3
    
    ctype=ctypes{ii};
    f=tester_crand(L,1);
    g=tester_crand(L,1);
    
    h1=pconv(f,g,ctype);
    h2=ref_pconv(f,g,ctype);
    
    res=norm(h1-h2);
    [test_failed,fail]=ltfatdiditfail(res,test_failed);        
    s=sprintf('PCONV %3s  L:%3i %0.5g %s',ctype,L,res,fail);
    disp(s);
    
  end;
end;


