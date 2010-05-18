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
  
% Parameters
sigLen = 450;
nbTimePosition = 11;
winType = 'hann'; % type of window
width = round(linspace(32,256,nbTimePosition)); % width of the windows
M=width;
a=round(width/4);

%timePosition = round(cumsum(width/4)); % position in time of the windows
timePosition = cumsum(a); % position in time of the windows
timePosition = timePosition - timePosition(1);

win = cell(nbTimePosition, 1);
for k=1:nbTimePosition
  win{k}=firwin(winType,width(k));
end



sig = rand(sigLen,1)-0.5; % random signal

% test using dual windows
wind = nsgabdual(win,a,sigLen);
[c,Ls] = nsdgt(sig,win,a,width);
sigd=insdgt(c,wind,a,sigLen);

res=norm(sigd-sig)/ norm(sig);
[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf(['DUAL  %0.5g %s\n'],res,fail);

% test using tight windows
wint = nsgabtight(win,a,sigLen);
[ct,Ls] = nsdgt(sig,wint,a,width);
sigt=insdgt(ct,wint,a,sigLen);

res=norm(sigt-sig)/ norm(sig);
[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf(['TIGHT %0.5g %s\n'],res,fail);

[A,B]=nsgabframebounds(win,a,sigLen);
