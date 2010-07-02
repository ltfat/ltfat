function test_failed=test_nsdgt()
%TEST_DGT Simple test of nsdgt and associated functions
%  Usage: testnsdgt()
% 
%  This function checks the exact reconstruction (up to numeric precision)
%  of the functions nsdgt and insdgt, when using dual windows computed with
%  nsgabdual, or tight windows computed with nsgabtight.
%
%  This test is done on a single short random signal, for only one given set
%  of windows.
%  A more systematic testing would be required for a complete validation of
%  these functions (in particular for inclusion of the functions in LTFAT)

%  Author: Florent Jaillet, 2009-05

test_failed=0;

disp(' ===============  TEST_NSDGT ================');
  
for ii=1:3
  % Parameters
  sigLen = 450;
  nbTimePosition = 11;
  winType = 'hann'; % type of window
  width = round(linspace(32,256,nbTimePosition)); % width of the windows
  M=width;
  if ii==2
    % Make M longer than the window length.
    M=M+2;
  end;
  
  if ii==3
    % Make M shorter than the window length.
    M=M-2;
  end;
  
  
  a=round(width/4);
    
  win = cell(nbTimePosition, 1);
  for k=1:nbTimePosition
    win{k}=firwin(winType,width(k));
  end



  sig = randn(sigLen,1)-0.5; % random signal
  
  for jj=1:2
    if jj==1
      winn = win;
      wind = nsgabdual(win,a,sigLen);
      wintype='DUAL ';
    else    
      wind = nsgabtight(win,a,sigLen);
      winn = wind;
      wintype='TIGHT';
    end;
    
    [c,Ls] = nsdgt(sig,winn,a,M);
    sigd=insdgt(c,wind,a,sigLen);
    
    res=norm(sigd-sig)/ norm(sig);
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    fprintf(['%s %0.5g %s\n'],wintype,res,fail);
    
  end;
  
end;

