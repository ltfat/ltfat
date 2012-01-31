function test_failed=test_pbspline

Lr=[15,16,18,20];
ar=[ 3, 4, 6, 5];
or=[1, 1.5, 2,3];

%btypes={'ed','xd','stard','ec','xc','starc'};
btypes={'ed','xd','stard'};
centtypes={'wp','hp'};

test_failed=0;

disp(' ===============  TEST_PBSPLINE ============');

for ii=1:length(Lr)
  L=Lr(ii);
  a=ar(ii);
  N=L/a;
  
  for jj=1:length(or)
    order=or(jj);
    
    for kk=1:numel(btypes)
      btype=btypes{kk};
      
      for ll=1:2
        centstring=centtypes{ll};
        
        [g,nlen]=pbspline(L,order,a,btype,centstring);
        
        A=zeros(L,1);
        
        for n=0:N-1
          A=A+circshift(g,n*a);
        end;
        
        res=max(abs(A-1/sqrt(a)));
        [test_failed,fail]=ltfatdiditfail(res,test_failed);        
        s=sprintf('PBSPLINE PU   %2s %s L:%3i a:%3i o:%3.5g %0.5g %s', ...
                  btype,centstring,L,a,order,res,fail);
        disp(s);
        
        gcutextend=middlepad(middlepad(g,nlen,centstring),L,centstring);
        
        res=norm(g-gcutextend);

        [test_failed,fail]=ltfatdiditfail(res,test_failed);        
        s=sprintf('PBSPLINE NLEN %2s %s L:%3i a:%3i o:%3.5g %0.5g %s', ...
                  btype,centstring,L,a,order,res,fail);
        disp(s);

        
      end;
      
    end;        
  end;
  
end;    

