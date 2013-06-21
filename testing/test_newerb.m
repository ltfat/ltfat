
testcase=2
isuniform=0;
switch testcase
  case 1
    [g,a]=erbfilters(16000,4000,'fractional');
    L=filterbanklength(4000,a);
    isreal=1;
  case 2
    [g,a]=erbfilters(16000);
    L=filterbanklength(4000,a);
    isreal=1;
  case 3
    [g,a]=erbfilters(16000,'uniform');
    L=filterbanklength(6000,a);
    isreal=1;
    isuniform=1;
end;

% Test it
if 0
    f=randn(L,1);   
else
    f=[1;zeros(L-1,1)];
    if 0
        ff=fft(f);
        ff(2000)=0;
        ff(2001)=0;
        ff(2002)=0;
        f=ifft(ff);
    end;
end;



% Inspect it: Dual windows, frame bounds and the response
disp('Frame bounds:')
[A,B]=filterbankrealbounds(g,a,L);
A
B
B/A
filterbankresponse(g,a,L,'real','plot');
gd=filterbankrealdual(g,a,L);

if isuniform
    c=ufilterbank(f,g,a);
else
    c=filterbank(f,g,a);
end;
r=ifilterbank(c,gd,a);
if isreal
    r=2*real(r);
end;

disp('Reconstruction:')
norm(f-r)


gt=filterbankrealtight(g,a,L);
if isuniform
    ct=ufilterbank(f,gt,a);
else
    ct=filterbank(f,gt,a);
end;
rt=ifilterbank(ct,gt,a);
if isreal
    rt=2*real(rt);
end;

disp('Reconstruction tight:')
norm(f-r)
