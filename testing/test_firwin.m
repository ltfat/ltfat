function test_failed=test_firwin
%TEST_FIRWIN  Test the firwin windows
%
%  This test script verifies the properties listed in the help of firwin
  
  
allwins = getfield(arg_firwin,'flags','wintype');

test_failed=0;

disp(' ===============  TEST_FIRWIN ================');

for L=[18,19,20,21]

  for cent=0:1
    if cent==0
      centtype='wp';
    else
      centtype='hp';
    end;
    for ii=1:length(allwins);
      winname=allwins{ii};
      
      [g,info]=firwin(winname,L,centtype);
      
      res = 1-isevenfunction(fir2long(g,2*L),centtype);
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      
      s=sprintf(['SYMM %10s %s L: %i %0.5g %s'],winname,centtype,L,res,fail);
      disp(s);
      
      if cent==0
        res=1-g(1);
        
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
        
        s=sprintf(['PEAK %10s %s L: %i %0.5g %s'],winname,centtype,L,res,fail);
        disp(s);
        
        
      end;
      
      if mod(L,2)==0 
        if info.ispu
          gpu=g+fftshift(g);
          res=norm(gpu-gpu(1)*ones(L,1));
          
          [test_failed,fail]=ltfatdiditfail(res,test_failed);
          
          s=sprintf(['PU   %10s %s L: %i %0.5g %s'],winname,centtype,L,res,fail);
          disp(s);
          
        end;
        
        if info.issqpu
          gpu=g.^2+fftshift(g.^2);
          res=norm(gpu-gpu(1)*ones(L,1));
          
          [test_failed,fail]=ltfatdiditfail(res,test_failed);
          
          s=sprintf(['SQPU %10s %s L: %i %0.5g %s'],winname,centtype,L,res,fail);
          disp(s);
          
        end;
      end;      
    end;
      
  end;
  
end;


