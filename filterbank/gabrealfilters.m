function [g,a,M2,cfreq,L] = gabrealfilters(a,M,Ls)

M2 = floor(M/2)+1;
cfreq = 2*(0:M2-1).'/M;

L = lcm(a,M);
L = ceil(Ls/L)*L;

g0 = fftshift(sqrt(L)*pgauss(L,1));
Lg = L;

g = cell(1,M2);
for kk = 0:M2-1
gtmp.H = g0;
gtmp.foff = kk*L/M-Lg/2;
gtmp.realonly = 0; 
gtmp.L = L;
g{kk+1} = gtmp;
%g3{kk+1}.H = comp_transferfunction(g2{kk+1},L);
%g3{kk+1}.foff = 0;
%g3{kk+1}.realonly = 0; 
end

%f = gspi;
%f = f(1:L);
%tic; [tgrad,fgrad,c_s]=filterbankphasegrad(f,g,a,L); PGtimeL = toc
%tic; sr=filterbankreassign(c_s,tgrad,fgrad,a,cfreq); RAtimeL = toc
%plotfilterbank(sr,a,'fc',22050*cfreq,'linabs');

%tic; [tgrad,fgrad,c_s]=filterbankphasegrad(f,g3,a,L); PGtimeL = toc
%tic; sr=filterbankreassign(c_s,tgrad,fgrad,a,cfreq); RAtimeL = toc
%subplot(313); plotfilterbank(sr,a,'fc',22050*cfreq,'linabs');

