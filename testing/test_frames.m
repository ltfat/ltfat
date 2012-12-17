function test_failed=test_frames
%TEST_FRAMES  Test the frames methods

test_failed=0;
  
disp(' ===============  TEST_FRAMES ================');


Fr=cell(1,26);

L=200;
Fr{1}  = frame('dgt','gauss',10,20);
Fr{2}  = frame('dgtreal','gauss',10,20);
Fr{3}  = frame('dwilt','gauss',20);
Fr{4}  = frame('wmdct','gauss',20);
Fr{5}  = frame('gen',crand(200,300),20);

Fr{6}  = frametight(frame('dgt','gauss',10,20));
Fr{7}  = frametight(frame('dgtreal','gauss',10,20));
Fr{8}  = frametight(frame('dwilt','gauss',20));
Fr{9}  = frametight(frame('wmdct','gauss',20));
Fr{10} = frametight(frame('gen',crand(200,300),20));

Fr{11} = frame('dft');
Fr{12} = frame('dcti');
Fr{13} = frame('dctii');
Fr{14} = frame('dctiii');
Fr{15} = frame('dctiv');
Fr{16} = frame('dsti');
Fr{17} = frame('dstii');
Fr{18} = frame('dstiii');
Fr{19} = frame('dstiv');

gfilt={randn(30,1),randn(20,1),randn(15,1),randn(10,1)};
Fr{20} = frame('ufilterbank',    gfilt,3,4);
Fr{21} = frame('ufilterbankreal',gfilt,3,4);

Fr{22} = frame('dgt','gauss',4,6,'lt',[1 2]);
Fr{23} = frame('identity');
Fr{24} = frame('fusion',[1 1],Fr{1},Fr{1});

%Fr{25} = frame('filterbank',     gfilt,[4 3 2 2],4);
%Fr{26} = frame('filterbankreal', gfilt,[4 3 2 2],4);



f=randn(L,1);

for ii=1:numel(Fr)
  F=Fr{ii};
  
  if isempty(F)
    continue;
  end;
  
  Fd=framedual(F);
  
  c=frana(F,f);
  r=frsyn(Fd,c);
  res=norm(r(1:L)-f);
  
  [test_failed,fail]=ltfatdiditfail(res,test_failed);
  s=sprintf(['FRAMES DUAL REC       frameno:%3i %s %0.5g %s'],ii,F.type,res,fail);    
  disp(s);  
  
  F2=frameaccel(F,L);
  F2d=frameaccel(Fd,L);

  c=frana(F2,f);
  r=frsyn(F2d,c);
  res=norm(r(1:L)-f);
  
  [test_failed,fail]=ltfatdiditfail(res,test_failed);
  s=sprintf(['FRAMES ACCEL DUAL REC frameno:%3i %s %0.5g %s'],ii,F.type,res,fail);    
  disp(s);

  [A,B]=framebounds(F,L);
  
  if ~any(strcmp(F.type,{'dgtreal','fftreal','ufilterbankreal','filterbankreal'}))
    LL=framelength(F,L);
    G=framematrix(F,LL);
    res=norm(c-G'*postpad(f,LL));
    
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    s=sprintf(['FRAMES ANA MATRIX     frameno:%3i %s %0.5g %s'],ii,F.type,res,fail);    
    disp(s);
    
    res=norm(frsyn(F,c)-G*c);
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    s=sprintf(['FRAMES SYN MATRIX     frameno:%3i %s %0.5g %s'],ii,F.type,res,fail);    
    disp(s);
    
  end;

end;
