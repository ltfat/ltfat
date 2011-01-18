function test_failed=test_firwin
%TEST_FIRWIN  Test the firwin windows
%
%  This test script verifies the properties listed in the help of firwin
  
  
PUwins = {'hann','halfsquare','tria'};

orthwins = {'sine','sqrthalfsquare','sqrttria','ogg'};

otherwins = {'square','hamming','blackman','nuttall'};

allwins = [PUwins,orthwins,otherwins];

test_failed=0;

disp(' ===============  TEST_FIRWIN ================');

for L=[18,19]

  for cent=0:1
    if cent==0
      centtype='wp';
    else
      centtype='hp';
    end;
    for ii=1:length(allwins);
      winname=allwins{ii};
      
      g=firwin(winname,L,centtype);
      
      res = 1-iseven(fir2long(g,2*L),centtype);
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      
      s=sprintf(['SYMM %10s %s L: %i %0.5g %s'],winname,centtype,L,res,fail);
      disp(s);
      
      if cent==0
        res=1-g(1);
        
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
        
        s=sprintf(['PEAK %10s %s L: %i %0.5g %s'],winname,centtype,L,res,fail);
        disp(s);
        
        
      end;
      
      
    end;
  
    % Only run this test for even length windows.
    if mod(L,2)==0    
      for ii=1:length(PUwins);
        winname=PUwins{ii};
        
        g=firwin(winname,L,centtype);
        
        gpu=g+fftshift(g);
        res=norm(gpu-gpu(1)*ones(L,1));
        
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
        
        s=sprintf(['PU   %10s %s L: %i %0.5g %s'],winname,centtype,L,res,fail);
        disp(s);
        
      end;
    end;
    
  end;
  
  % Only run this test for even length windows.
  if mod(L,2)==0
    for ii=1:length(orthwins);
      winname=orthwins{ii};
      
      g=firwin(winname,L,'2');
      
      [A,B] = wilbounds(g,L/2);
      
      res=abs(A-1)+abs(B-1);
      
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      
      s=sprintf(['ORTH %10s L: %i A:%0.5g B:%0.5g %0.5g %s'],winname,L,A,B,res,fail);
      disp(s);
      
    end;
  end;
  
end;
