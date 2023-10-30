function test_failed = test_hermbasis()

test_failed = 0;

Lr=[100,101,102,103];

test_failed=0;

highordercoef=0.80;

disp(' ===============  TEST_HERMBASIS ==========');

for ii=1:length(Lr)
  
    for n=2:4
        
        L=Lr(ii);

        highorder=round(L*highordercoef);

        [H,D]=hermbasis(L,n);
        
        r1=(H*H')-eye(L);
        res=norm(r1,'fro');
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
        
        s=sprintf('HERMBASIS orth L:%3i n:%3i %0.5g %s',L,n,res,fail);
        disp(s);
        
        f=tester_crand(L,1);
               
        res=norm(dft(f)-H*diag(D)*H'*f);
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
        
        s=sprintf('HERMBASIS DFT  L:%3i n:%3i %0.5g %s',L,n,res,fail);
        disp(s);
        
        [H,D]=pherm(L,0:highorder-1,'fast');

        res=norm(dft(H)-H*diag(D));
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
        
    end;
        
    s=sprintf('PHERM     DFT  L:%3i %0.5g %s',L,res,fail);
    disp(s);
    
    [H,D]=pherm(L,0:highorder-1,'fast','qr');
    
    res=norm(dft(H)-H*diag(D));
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    
    s=sprintf('PHERM QR  DFT  L:%3i %0.5g %s',L,res,fail);
    disp(s);
    
    r1=(H'*H)-eye(highorder);
    res=norm(r1,'fro');
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    
    s=sprintf('PHERM QR orth  L:%3i %0.5g %s',L,res,fail);
    disp(s);
    
    
    [H,D]=pherm(L,0:highorder-1,'fast','polar');
    
    res=norm(dft(H)-H*diag(D));
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    
    s=sprintf('PHERM POL DFT  L:%3i %0.5g %s',L,res,fail);
    disp(s);
    
    r1=(H'*H)-eye(highorder);
    res=norm(r1,'fro');
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    
    s=sprintf('PHERM POL orth L:%3i %0.5g %s',L,res,fail);
    disp(s);
            
end;

