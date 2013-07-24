function test_failed=test_frames
%TEST_FRAMES  Test the frames methods

test_failed=0;
  
disp(' ===============  TEST_FRAMES ================');
global LTFAT_TEST_TYPE;

tolchooser.double=1e-9;
tolchooser.single=2e-4;
tolerance = tolchooser.(LTFAT_TEST_TYPE);

% Iterative algorithms need a bigger tolerance
tolchooseriter.double=1e-1;
tolchooseriter.single=1e-1;
toleranceiter = tolchooseriter.(LTFAT_TEST_TYPE);

Fr=cell(1,32);

L=200;
Fr{1}  = frame('dgt','gauss',10,20);
Fr{2}  = frame('dgtreal','gauss',10,20);
Fr{3}  = frame('dwilt','gauss',20);
Fr{4}  = frame('wmdct','gauss',20);
Fr{5}  = frame('gen',tester_crand(200,300),20);

Fr{6}  = frametight(frame('dgt','gauss',10,20));
Fr{7}  = frametight(frame('dgtreal','gauss',10,20));
Fr{8}  = frametight(frame('dwilt','gauss',20));
Fr{9}  = frametight(frame('wmdct','gauss',20));
Fr{10} = frametight(frame('gen',tester_crand(200,300),20));

Fr{11} = frame('dft');
Fr{12} = frame('dcti');
Fr{13} = frame('dctii');
Fr{14} = frame('dctiii');
Fr{15} = frame('dctiv');
Fr{16} = frame('dsti');
Fr{17} = frame('dstii');
Fr{18} = frame('dstiii');
Fr{19} = frame('dstiv');

gfilt={tester_rand(30,1),tester_rand(20,1),tester_rand(15,1),tester_rand(10,1)};
Fr{20} = frame('ufilterbank',    gfilt,3,4);
Fr{21} = frame('ufilterbankreal',gfilt,3,4);

Fr{22} = frame('dgt','gauss',4,6,'lt',[1 2]);
Fr{23} = frame('identity');
Fr{24} = frame('fusion',[1 1],Fr{1},Fr{1});
Fr{25} = frametight(frame('dgt','hamming',10,20));
Fr{26} = frametight(frame('wmdct','hamming',20));

g={randn(30,1),randn(50,1),randn(70,1),randn(90,1)};
a=[20,40,60,80];
M=[30,50,70,100];

Fr{27} = frametight(frame('nsdgt',g,a,M));
Fr{28} = frametight(frame('unsdgt',g,a,100));
Fr{29} = frametight(frame('nsdgtreal',g,a,M));
Fr{30} = frametight(frame('unsdgtreal',g,a,100));

Fr{31} = frametight(frame('dftreal'));
Fr{32} = frame('fwt','db4',5);

% The tensor frame implementation is currenly broken
%Fr{33} = frame('tensor',Fr{11});


%Fr{31} = frame('filterbank',     gfilt,[4 3 2 2],4);
%Fr{32} = frame('filterbankreal', gfilt,[4 3 2 2],4);



f=tester_rand(L,1);

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

  % Test that framebounds are able to run, not actual resting is done on
  % the values.
  [A,B]=framebounds(F,L);
  
  %% Test iterative analysis and synthesis
  r=frsyniter(F,c);
  
  res=norm(r(1:L)-f);
  [test_failed,fail]=ltfatdiditfail(res,test_failed,toleranceiter);
  s=sprintf(['FRSYNITER             frameno:%3i %s %0.5g %s'],ii,F.type,res,fail);    
  disp(s);
  
  c2=franaiter(Fd,f);
  res=norm(c2-c);
  [test_failed,fail]=ltfatdiditfail(res,test_failed,toleranceiter);
  s=sprintf(['FRANAITER             frameno:%3i %s %0.5g %s'],ii,F.type,res,fail);    
  disp(s);  
  
  %% Test matrix representations
  if (~F.realinput)
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
  
  %% Test the frame multipliers: test framemul, framemuladj and
  %% iframemul
  if F.realinput
      m=1+0.01*tester_rand(size(c,1),1);
  else
      m=1+1i+0.01*tester_crand(size(c,1),1);
  end;
  ff=framemul(f,F,Fd,m);
  fr=iframemul(ff,F,Fd,m);
  res=norm(f-fr(1:L))/norm(f);
  [test_failed,fail]=ltfatdiditfail(res,test_failed,tolerance);
  s=sprintf('IFRAMEMUL             frameno:%3i %s %0.5g %s',ii, ...
            F.type,res,fail);
  disp(s);
  
end;
