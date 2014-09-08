function test_failed=test_ptpfun()
test_failed = 0;

% First, just test if functions run with various input parameters

M = 10;
a = 5;
L = 20;
incrange = 0:4:80;

g = ptpfun(L,[1,-1]);
g = ptpfun(L,[1,-1,9]);
g = ptpfun(L,[1,-1,9],'inf');

gd = ptpfundual(L,[1,-1],a,M);
gd = ptpfundual(L,[1,-1],a,M,10);

[gd,nlen] = ptpfundual(L,[1,-1],a,M,'inf');

[gd,nlen] = ptpfundual(L,[1,-1],a,M,'inf','matchscale');


% This should fail, but be caught by some of the input checks.
% We will test if the error message starts with function name in allcaps
% followed by a colon e.g. PTPFUN:

% Too short w
try
    g = ptpfun(L,[1]);
    gd = ptpfundual(L,[1],a,M);
    % We should have failed in g
    test_failed = test_failed + 1;
    failstr = 'FAILED';
catch
    err = lasterror;
    [test_failed,failstr]=dititfailedcorrectly(err.message,'ptpfun',test_failed);
end
fprintf('PTPFUN Too short w test %s\n',failstr);

% Only pos weights
try
    g = ptpfun(L,[1,1]);
    gd = ptpfundual(L,[1,1],a,M);
    % We should have failed in g
    test_failed = test_failed + 1;
    failstr = 'FAILED';
catch
    err = lasterror;
    [test_failed,failstr]=dititfailedcorrectly(err.message,'ptpfun',test_failed);
end
fprintf('PTPFUN Only positive w test %s\n',failstr);

% Only neg weights
try
    g = ptpfun(L,[-1,-1]);
    gd = ptpfundual(L,[-1,-1],a,M);
    % We should have failed in g
    test_failed = test_failed + 1;
    failstr = 'FAILED';
catch
    err = lasterror;
    [test_failed,failstr]=dititfailedcorrectly(err.message,'ptpfun',test_failed);
end
fprintf('PTPFUN Only negative w test %s\n',failstr);

% One zero in weights
try
    g = ptpfun(L,[-1,0,1]);
    gd = ptpfundual(L,[-1,0,1],a,M);
    % We should have failed in g
    test_failed = test_failed + 1;
    failstr = 'FAILED';
catch
    err = lasterror;
    [test_failed,failstr]=dititfailedcorrectly(err.message,'ptpfun',test_failed);
end
fprintf('PTPFUN Zero in w test %s\n',failstr);


% Test if ptpfun and ptpfundual indeed fulfill the Waxler-Raz conditions
% (are dual windoes up to scaling) for a range of inc parameter
wcell = {[-1,1],[1,-1,3,4],[7,8,-3]};

for w = wcell
for inc =incrange
   g = ptpfun(L,w{1});
   gd = ptpfundual(L,w{1},a,M);
   [~,err] = gabdualnorm(g,gd,a,M,L);
   [test_failed,fail]=ltfatdiditfail(err,test_failed);
   fprintf('PTPFUN IS DUAL L=%i,a=%i,M=%i, inc=%i %s\n',L,a,M,inc,fail);
end
end

% Test ptpfun and ptpfundual individually (each uses canonical dual window)
f = tester_crand(L,1);

for w = wcell
   g = ptpfun(L,w{1});
   c = dgt(f,g,a,M);
   fhat = idgt(c,{'dual',g},a);
   res = norm(f-fhat);
   [test_failed,fail]=ltfatdiditfail(res,test_failed);
   fprintf('PTPFUN REC L=%i,a=%i,M=%i, %s\n',L,a,M,fail);
end

for w = wcell
   g = ptpfundual(L,w{1},a,M);
   c = dgt(f,g,a,M);
   fhat = idgt(c,{'dual',g},a);
   res = norm(f-fhat);
   [test_failed,fail]=ltfatdiditfail(res,test_failed);
   fprintf('PTPFUNDUAL REC L=%i,a=%i,M=%i, %s\n',L,a,M,fail);
end

% Test ptpfun and ptpfundual properly scaled
for w = wcell
  for inc =incrange
   g = ptpfun(L,w{1});
   gd = ptpfundual(L,w{1},a,M,inc,'matchscale');
   c = dgt(f,g,a,M);
   fhat = idgt(c,gd,a);
   res = norm(f-fhat);
   [test_failed,fail]=ltfatdiditfail(res,test_failed);
   fprintf('PTPFUN PTPFUNDUAL REC L=%i,a=%i,M=%i, inc=%i %s\n',L,a,M,inc,fail);
  end
end



function [test_failed,failstr]=dititfailedcorrectly(errmsg,fname,test_failed)

if isempty(strfind(errmsg,strcat(upper(fname),':')))
    test_failed = test_failed + 1;
    failstr = 'FAILED';
else
    failstr = '';
end

