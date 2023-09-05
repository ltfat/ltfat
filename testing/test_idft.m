function test_failed=test_idft
Lr=[1, 19, 20];


test_failed=0;

disp(' ===============  TEST_IDFT ==============');

for jj=1:length(Lr)
  L=Lr(jj);
    for n = 1:2
    
    if (n==1)
       type = 'complex';
       c=tester_crand(L,1);
    elseif (n==2)
       type = 'real';
       c=tester_rand(L,1);      
    end
    
    f1=idft(c);
    f2=ref_idft(c);
    
    res=norm(f1-f2);
    [test_failed,fail]=ltfatdiditfail(res,test_failed);        
    s=sprintf('IDFT %6s  L:%3i %0.5g %s',type,L,res,fail);
    disp(s);
    end
  end;
%end;


