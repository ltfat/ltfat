function test_failed=test_dft
Lr=[1, 19, 20];


test_failed=0;

disp(' ===============  TEST_DFT ==============');

for jj=1:length(Lr)
  L=Lr(jj);
    for n = 1:2
    
    if (n==1)
       type = 'complex';
       f=tester_crand(L,1);
    elseif (n==2)
       type = 'real';
       f=tester_rand(L,1);      
    end
    
    c1=dft(f);
    c2=ref_dft(f);
    
    res=norm(c1-c2);
    [test_failed,fail]=ltfatdiditfail(res,test_failed);        
    s=sprintf('DFT %6s  L:%3i %0.5g %s',type,L,res,fail);
    disp(s);
    end
  end;
%end;


