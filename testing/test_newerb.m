
if 0
    [g,a]=erbfilters(16000,4000,'fractional');
    L=filterbanklength(4000,a);
end;

if 0
    [g,a]=erbfilters(16000);
    L=filterbanklength(4000,a);
end;

if 1
    [g,a]=erbfilters(16000,'uniform');
    L=filterbanklength(6000,a);
end;



% Inspect it: Dual windows, frame bounds and the response
disp('Frame bounds:')
[A,B]=filterbankrealbounds(g,a,L);
A
B
B/A
filterbankresponse(g,a,L,'real','plot');
gd=filterbankrealdual(g,a,L);

% Test it
if 0
    f=randn(L,1);   
else
    f=[1;zeros(L-1,1)];
end;
c=filterbank(f,g,a);
rc=ifilterbank(c,gd,a);
r=2*real(rc);
disp('Reconstruction:')
norm(f-r)

