function test_failed=test_gabfirtight
%TEST_GABFIRDUAL  Test GABFIRTIGHT

      
Lr=[24,16,144,144];
ar=[ 4, 4,  8, 12];
Mr=[ 12, 16, 72, 48];

test_failed=0;
tolerance=1e-5;

disp(' ===============  TEST_GABFIRTIGHT ================');

disp('--- Used subroutines ---');


addpath([ltfatbasepath, 'thirdparty', filesep, 'unlocbox']);
if exist('init_unlocbox.m', 'file')
    init_unlocbox;
else
    disp('unlocbox not found. please initialize git submodule via git submodule update --init.')
end

for ii=1:length(Lr);

  L=Lr(ii);
  
  M=Mr(ii);
  a=ar(ii);
  b=L/M;
  N=L/a;
  c=gcd(a,M);
  d=gcd(b,N);
  p=a/c;
  q=M/c;
  

  for rtype=1:2
      
    if rtype==1
      rname='REAL ';	
      g=tester_rand(L,1);
    else
      rname='CMPLX';	
      g=tester_crand(L,1);
    end;
 
    global LTFAT_TEST_TYPE;
    if strcmpi(LTFAT_TEST_TYPE,'single')
        C = gabframebounds(g,a,M);
        while C>1e3
%             warning(sprintf(['The frame is too badly conditioned '...
%                              'for single precision. Cond. num. %d. '...
%                              ' Trying again.'],C));
                         
                         if rtype==1
                             rname='REAL ';
                             g=tester_rand(L,1);
                         else
                             rname='CMPLX';
                             g=tester_crand(L,1);
                         end;
                         C = gabframebounds(g,a,M);
        end
    end
    
    gt=gabfirtight(L, g,a,M);


    for W=1:3
          
      if rtype==1
        f=tester_rand(L,W);
      else
        f=tester_crand(L,W);
      end;

     
      
      
      % --- Test reconstruction of IDGT using the tight window. ---
      
      res=norm(f-idgt(dgt(f,gt,a,M),gt,a),'fro');
      [test_failed,fail]=ltfatdiditfail(res,test_failed, tolerance);
      s=sprintf(['TIG %s L:%3i W:%2i a:%3i b:%3i c:%3i d:%3i p:%3i q:%3i ' ...
                 '%0.5g %s'],rname,L,W,a,b,c,d,p,q,res,fail);
      disp(s);
      

    end;

  end;  

end;

close_unlocbox;
rmpath(([ltfatbasepath, 'thirdparty', filesep, 'unlocbox']));
