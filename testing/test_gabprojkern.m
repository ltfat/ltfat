function test_failed = test_gabprojkern
test_failed = 0;
tol = 1e-12;

test = gspi;
L = 22050;
x = test(1:L);

g = 'gauss';
a = 16;
M = 1024;

% build frame analysis coefficients
c = dgt(x,g,a,M);

% test if this is indeed a projection
cout = gabprojkern(c,g,a);
cout2 = gabprojkern(cout,g,a);

res = norm(cout-cout2);

[test_failed,fail]=ltfatdiditfail(res,test_failed,tol);
fprintf('GABPROJKERN COMPLEX %s\n',fail);

% build frame analysis coefficients
c = dgtreal(x,g,a,M);

% test if this is indeed a projection
cout = gabprojkern(c,g,a,M,'real');
cout2 = gabprojkern(cout,g,a,M,'real');

res = norm(cout-cout2);

[test_failed,fail]=ltfatdiditfail(res,test_failed,tol);
fprintf('GABPROJKERN REAL %s\n',fail);