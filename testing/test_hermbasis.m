Lr=[9,10,12];

test_failed=0;

disp(' ===============  TEST_HERMBASIS ==========');

for ii=1:length(Lr)
  
    for n=2:4
        
        L=Lr(ii);
        
        [H,D]=hermbasis(L,n);
        
        r1=(H*H')-eye(L);
        res=norm(r1,'fro');
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
        
        s=sprintf('HERMBASIS orth L:%3i n:%3i %0.5g %s',L,n,res,fail);
        disp(s);
        
        f=crand(L,1);
               
        res=norm(dft(f)-H*diag(D)*H'*f);
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
        
        s=sprintf('HERMBASIS DFT  L:%3i n:%3i %0.5g %s',L,n,res,fail);
        disp(s);

        [H,D]=pherm(L,0:round(L*.8)-1);

        res=norm(dft(H)-H*diag(D));
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
        
        s=sprintf('PHERM     DFT  L:%3i n:%3i %0.5g %s',L,n,res,fail);
        disp(s);

        
    end;
        
end;

