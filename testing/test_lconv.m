function test_failed=test_lconv
Lr=[1,19,20];

ctypes={'default','r','rr'};

test_failed=0;

disp(' ===============  TEST_LCONV ==============');

for jj=1:length(Lr)
  L=Lr(jj);
  for ii=1:3
    for type = {'real','complex'}
    ctype=ctypes{ii};
    
    if strcmp(type{1},'complex')
       f=tester_crand(L,1);
       g=tester_crand(L,1);
    else
       f=tester_rand(L,1);
       g=tester_rand(L,1);        
    end
    
    h1=lconv(f,g,ctype);
    h2=ref_lconv(f,g,ctype);
    
    res=norm(h1-h2);
    [test_failed,fail]=ltfatdiditfail(res,test_failed);        
    s=sprintf('LCONV %3s %6s  L:%3i %0.5g %s',ctype,type{1},L,res,fail);
    disp(s);
    end
  end;
end;


