function test_failed=test_frames
%TEST_FRAMES  Test the frames methods

test_failed=0;
  
disp(' ===============  TEST_FRAMES ================');


Fr=cell(1,2);

L=200;
Fr{1} =newframe('dgt','gauss','dual',10,20);
Fr{2} =newframe('dgtreal','gauss','dual',10,20);
Fr{3} =newframe('dwilt','gauss','dual',20);
Fr{4} =newframe('wmdct','gauss','dual',20);
Fr{5} =newframe('gen',crand(200,300),'dual',20);

Fr{6} =newframe('dgt','gauss','tight',10,20);
Fr{7} =newframe('dgtreal','gauss','tight',10,20);
Fr{8} =newframe('dwilt','gauss','tight',20);
Fr{9} =newframe('wmdct','gauss','tight',20);
Fr{10} =newframe('gen',crand(200,300),'tight',20);

Fr{11} =newframe('dft');
Fr{12} =newframe('fft');
Fr{13} =newframe('fftreal',200);
Fr{14} =newframe('dcti');
Fr{15} =newframe('dctii');
Fr{16} =newframe('dctiii');
Fr{17}=newframe('dctiv');
Fr{18}=newframe('dsti');
Fr{19}=newframe('dstii');
Fr{20}=newframe('dstiii');
Fr{21}=newframe('dstiv');

gfilt={randn(30,1),randn(20,1),randn(15,1),randn(10,1)};
Fr{22}=newframe('ufilterbank',    gfilt,'dual',3,4);
Fr{23}=newframe('ufilterbankreal',gfilt,'dual',3,4);

Fr{24}=newframe('dgt','gauss','dual',4,6,'lt',[1 2]);
Fr{25}=newframe('identity');

%Fr{19}=newframe('filterbank',     gfilt,[],[4 3 2 2],4);
%Fr{20}=newframe('filterbankreal', gfilt,[],[4 3 2 2],4);



f=randn(L,1);

for ii=1:numel(Fr)
  F=Fr{ii};
  
  if isempty(F)
    continue;
  end;
  
  FT=frametype(F);
  
  c=frana(F,f);
  r=frsyn(F,c);
  res=norm(r(1:L)-f);
  
  [test_failed,fail]=ltfatdiditfail(res,test_failed);
  s=sprintf(['FRAMES RECONSTRUCTION %s %0.5g %s'],FT,res,fail);    
  disp(s);  

  F2=frameaccel(F,L);

  c=frana(F2,f);
  r=frsyn(F2,c);
  res=norm(r(1:L)-f);
  
  [test_failed,fail]=ltfatdiditfail(res,test_failed);
  s=sprintf(['FRAMES ACCEL RECON    %s %0.5g %s'],FT,res,fail);    
  disp(s);

  cadj=frsynadj(F,f);
  radj=franaadj(F,cadj);
  res=norm(radj(1:L)-f);
  
  [test_failed,fail]=ltfatdiditfail(res,test_failed);
  s=sprintf(['FRAMES ADJ RECON      %s %0.5g %s'],FT,res,fail);    
  disp(s);  

  
  if ~any(strcmp(FT,{'dgtreal','fftreal','ufilterbankreal','filterbankreal'}))
    LL=framelength(F,L);
    G=franamat(F,LL);
    res=norm(c-G'*postpad(f,LL));
    
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    s=sprintf(['FRAMES MATRIX         %s %0.5g %s'],FT,res,fail);    
    disp(s);
    
    res=norm(franaadj(F,c)-G*c);
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    s=sprintf(['FRAMES ANA ADJOINT    %s %0.5g %s'],FT,res,fail);    
    disp(s);
    
  end;

end;
