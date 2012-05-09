function test_failed=test_nsdgt()
%TEST_NSDGT Simple test of nsdgt and associated functions
%  Usage: test_nsdgt()
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

a=[20,30,40];
M=[30,40,50];
L=sum(a);

f=randn(L,1);

g=cell(3,1);
for ii=1:3
  g{ii}=randn(M(ii),1);
end;

% ----- non-uniform dual and inversion -----

gd=nsgabdual(g,a,M);

c=nsdgt(f,g,a,M);
r=insdgt(c,gd,a);

res=norm(f-r);

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf(['NSDGT DUAL  %0.5g %s\n'],res,fail);

% ----- non-uniform inversion, real -----

c=nsdgtreal(f,g,a,M);
r=insdgtreal(c,gd,a,M);

res=norm(f-r);

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf(['NSDGTREAL DUAL  %0.5g %s\n'],res,fail);



% ----- tight and inversion -----------------
gt=nsgabtight(g,a,M);

ct=nsdgt(f,gt,a,M);
rt=insdgt(ct,gt,a);

res=norm(f-rt);

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf(['NSDGT TIGHT %0.5g %s\n'],res,fail);



a=[20,30,40];
M=50;
L=sum(a);

f=randn(L,1);

g=cell(3,1);
for ii=1:3
  g{ii}=randn(M,1);
end;

% ----- non-uniform dual and inversion -----

gd=nsgabdual(g,a,M);

c=unsdgt(f,g,a,M);
r=insdgt(c,gd,a);

res=norm(f-r);

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf(['UNSDGT DUAL  %0.5g %s\n'],res,fail);

% ----- non-uniform inversion, real -----

c=unsdgtreal(f,g,a,M);
r=insdgtreal(c,gd,a,M);

res=norm(f-r);

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf(['UNSDGTREAL DUAL  %0.5g %s\n'],res,fail);



% ----- tight and inversion -----------------
gt=nsgabtight(g,a,M);

ct=unsdgt(f,gt,a,M);
rt=insdgt(ct,gt,a);

res=norm(f-rt);

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf(['UNSDGT TIGHT %0.5g %s\n'],res,fail);

% ------ Reference DGT ----------------------

a1=3;
M1=4;
N=8;
a=a1*ones(1,N);
M=M1*ones(1,N);
L=a1*N;
f=crand(L,1);
g1=crand(L,1);
for ii=1:N
  g{ii}=g1;
end;

c     = nsdgt(f,g,a,M);
c_ref = dgt(f,g1,a1,M1,'timeinv');

res=norm(reshape(cell2mat(c),M1,N)-c_ref,'fro');

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf(['NSDGT REF %0.5g %s\n'],res,fail);

