function test_failed=test_frames
%TEST_FRAMES  Test the frames methods

test_failed=0;
  
disp(' ===============  TEST_DGT ================');


Fr=cell(1,2);

L=200;
Fr{1} =newframe('dgt','gauss','dual',10,20);
Fr{2} =newframe('dgtreal','gauss','dual',10,20);
Fr{3} =newframe('dwilt','gauss','dual',20);
Fr{4} =newframe('wmdct','gauss','dual',20);
Fr{5} =newframe('gen',crand(200,300),'dual',20);
Fr{6} =newframe('dft');
Fr{6} =newframe('fft');
Fr{6} =newframe('fftreal',200);
Fr{7} =newframe('dcti');
Fr{8} =newframe('dctii');
Fr{9} =newframe('dctiii');
Fr{10}=newframe('dctiv');
Fr{11}=newframe('dsti');
Fr{12}=newframe('dstii');
Fr{13}=newframe('dstiii');
Fr{14}=newframe('dstiv');

% The window support is not yeat good enough to handle filterbanks
%Fr{15}=newframe('ufilterbank',{randn(L,1),randn(L,1),randn(L,1),randn(L,1)},'dual',2);

f=randn(L,1);

for ii=1:numel(Fr)
  F=Fr{ii};
  
  if isempty(F)
    continue;
  end;
  
  FT=frametype(F);
  
  c=frana(F,f);
  r=frsyn(F,c);
  res=norm(r-f);
  
  [test_failed,fail]=ltfatdiditfail(res,test_failed);
  s=sprintf(['FRAMES RECONSTRUCTION %s %0.5g %s'],FT,res,fail);    
  disp(s);  

  F2=frameaccel(F,L);

  c=frana(F2,f);
  r=frsyn(F2,c);
  res=norm(r-f);
  
  [test_failed,fail]=ltfatdiditfail(res,test_failed);
  s=sprintf(['FRAMES ACCEL RECON    %s %0.5g %s'],FT,res,fail);    
  disp(s);

  if ~any(strcmp(FT,{'dgtreal','fftreal'}))
    G=framematrix(F,L);
    res=norm(c-G'*f);
    
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    s=sprintf(['FRAMES MATRIX         %s %0.5g %s'],FT,res,fail);    
    disp(s);
    
    res=norm(franaadj(F,c)-G*c);
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    s=sprintf(['FRAMES ANA ADJOINT    %s %0.5g %s'],FT,res,fail);    
    disp(s);
    
  end;

end;
