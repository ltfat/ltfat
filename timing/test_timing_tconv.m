%TEST_TIMING_DGT_FAC  Test timing factorization DGTs
%
%   This script test the timing SPREADADJs by comparing the results to
%   spreadadj in the main toolbox. Therefore, the correctness of
%   spreadadj must be verified first.


routinemax=2;

spfraction=.1;

test_failed=0;

disp('--- Used subroutines ---');


for rtype=1:2
  
  if rtype==1
    rname='REAL ';	
  else
    rname='CMPLX';	
  end;
  
  for sptype=1:2
    
    if sptype==1
      spname='FULL  ';	
    else
      spname='SPARSE';	
    end;
    
    for L=12:13

      if rtype==1
        if sptype==1
          coef1=rand(L,L);
          coef2=rand(L,L);
        else
          coef1=sprand(L,L,spfraction);
          coef2=sprand(L,L,spfraction);
        end;
      else
        if sptype==1
          coef1=crand(L,L);
          coef2=crand(L,L);
        else
          coef1=spcrand(L,L,spfraction);
          coef2=spcrand(L,L,spfraction);
        end;
      end;      
      
      ctwist=tconv(coef1,coef2);
      
      for rout=1:routinemax                        
        
        ctwist2=feval(['ref_tconv_',num2str(rout)],coef1,coef2);
                  
        rdiff=ctwist-ctwist2;
        
        res=norm(rdiff(:));      
        
        fail='';
        if res>10e-10
          fail='FAILED';
          test_failed=test_failed+1;
        end;
        
        s=sprintf('TWI %s %s %i L:%3i %0.5g %s',rname,spname,rout,L,res,fail);
        disp(s)
      end;

      
    end;

  end;  
  
end;


test_failed


