function test_failed=test_filterbank
%TEST_FILTERBANK test the filterbank codes
%  Usage: test_filterbank()
%
%  This function checks the exact reconstruction (up to numeric precision)
%  of the functions ufilterbank and ifilterbank, when using dual windows computed with
%  filterbankdual / filterbankrealdual, or tight windows computed with
%  filterbanktight / filterbankrealtight

test_failed=0;

disp(' ===============  TEST_FILTERBANK ================');

which comp_ufilterbank_fft

M=6;
a=3;
N=4;

L=a*N;

g=cell(1,M);
for ii=1:M
  g{ii}=crand(L,1);
end;

gd = filterbankdual(g,a);

f=crand(L,1);

c_u      = ufilterbank(f,g,a);
c_u_ref  = ref_ufilterbank(f,g,a);
c_nu     = filterbank(f,g,a);

%% check that filterbank and ufilterbank produce the same results.
res=0;
for m=1:M
  res=res+norm(c_nu{m}-c_u(:,m));  
end;

[test_failed,fail]=ltfatdiditfail(res,test_failed);
s=sprintf(['FB MATCH  %0.5g %s'],res,fail);
disp(s)

%% check that ufilterbank match its reference
res=norm(c_u-c_u_ref,'fro');

[test_failed,fail]=ltfatdiditfail(res,test_failed);
s=sprintf(['FB RES    %0.5g %s'],res,fail);
disp(s)


%% Check that ufilterbank is invertible using dual window
r=ifilterbank(c_u,gd,a);

res=norm(f-r);

[test_failed,fail]=ltfatdiditfail(res,test_failed);
s=sprintf(['FB DUAL   %0.5g %s'],res,fail);
disp(s)


%% Check that filterbanktight gives a tight filterbank
gt = filterbanktight(g,a);

c_ut = ufilterbank(f,gt,a);
r=ifilterbank(c_ut,gt,a);

res=norm(f-r);

[test_failed,fail]=ltfatdiditfail(res,test_failed);
s=sprintf(['FB TIGHT  %0.5g %s'],res,fail);
disp(s)

%% Check that filterbankbounds detect the tight frame
[AF,BF]=filterbankbounds(gt,a);

[test_failed,fail]=ltfatdiditfail(res,test_failed);
s=sprintf(['FB FB B   %0.5g %s'],BF-1,fail);
disp(s)

[test_failed,fail]=ltfatdiditfail(res,test_failed);
s=sprintf(['FB FB A   %0.5g %s'],AF-1,fail);
disp(s)

%% Check the real valued systems, dual

fr=rand(L,1);

gdreal=filterbankrealdual(g,a);

c_ur=ufilterbank(fr,g,a);
rreal=2*real(ifilterbank(c_ur,gdreal,a));

res=norm(fr-rreal);
[test_failed,fail]=ltfatdiditfail(res,test_failed);
s=sprintf(['FB RDUAL  %0.5g %s'],res,fail);
disp(s)


%% Check the real valued systems, tight

gtreal=filterbankrealtight(g,a);

ct     = ufilterbank(fr,gtreal,a);
rrealt = 2*real(ifilterbank(ct,gtreal,a));

res=norm(fr-rrealt);
[test_failed,fail]=ltfatdiditfail(res,test_failed);
s=sprintf(['FB RTIGHT %0.5g %s'],res,fail);
disp(s)

%% check filterbankrealbounds

[AF,BF]=filterbankrealbounds(gtreal,a);

[test_failed,fail]=ltfatdiditfail(res,test_failed);
s=sprintf(['FB FBR B  %0.5g %s'],BF-1,fail);
disp(s)

[test_failed,fail]=ltfatdiditfail(res,test_failed);
s=sprintf(['FB FBR A  %0.5g %s'],AF-1,fail);
disp(s)

