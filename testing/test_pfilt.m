function test_failed=test_pfilt
Lr =[27,100,213];

gr{1}=randn(20,1);
gr{2}=randn(21,1);
gr{3}=firfilter('hanning',19);
gr{4}=firfilter('hanning',20);
gr{5}=randn(4,1);
gr{6}=firfilter('hanning',20,'causal');
gr{7}=firfilter('hanning',20,'delay',13); 
gr{8}=firfilter('hanning',20,'delay',-13); 
gr{7}=firfilter('hanning',20,'delay',14); 
gr{8}=firfilter('hanning',20,'delay',-14); 
gr{9}=firfilter('hamming',19,.3);   
gr{10}=firfilter('hamming',19,.3,'real');
gr{11}=blfilter('hanning',.19);
gr{12}=blfilter('hanning',.2);
gr{13}=blfilter('hanning',.132304);
gr{14}=blfilter('hanning',.23,'delay',13);
gr{15}=blfilter('hamming',.23,.3);
gr{16}=blfilter('hamming',.23,.3,'real');
gr{17}=blfilter('hanning',2);

test_failed=0;

disp(' ===============  TEST_PFILT ==============');

disp('--- Used subroutines ---');

for ii=1:numel(gr)
  g=gr{ii};


  for a=1:3

      for jj=1:length(Lr)
          L=ceil(Lr(jj)/a)*a;
      
          for W=1:3
              
              for rtype=1:2
                  if rtype==1
                      rname='REAL ';	
                      f=tester_rand(L,W);
                  else
                      rname='CMPLX';	
                      f=tester_crand(L,W);
                  end;
                  
                                   
                  h2=ref_pfilt(f,g,a);
                  h1=pfilt(f,g,a);
				  
                   res=norm(h1-h2);
                  [test_failed,fail]=ltfatdiditfail(res,test_failed);        
                  s=sprintf('PFILT %3s  filtno:%3i L:%3i W:%3i a:%3i %0.5g %s',rname,ii,L,W,a,res,fail);
                  disp(s);
              end;
          end;
      end;
      
  end;
  
end;

