Ls = 44100;
[f,fs] = gspi;
f = f(1:Ls);

% Create an UNIFORM filterbank
a = 32;
M = 1024;

%% Regular uniform Gabor filters

[g,a,M,cfreq,L] = gabrealfilters(a,M,Ls,1);
%[g,a,fc] = erbfilters(fs,Ls,'uniform','spacing',1/10,'bwmul',1);

% We will not do no subsampling. 
% This is the main requirement for synchrosqueezing to work.
%a = 1;

[tgradtmp,fgradtmp,~,ctmp] = filterbankphasegrad(f,g,a,L);
 
M = size(ctmp,1);
N = length(ctmp{1});
c=zeros(M,N);
tgrad = c; fgrad = c;
for m=1:M    
    c(m,:)=ctmp{m};    
    tgrad(m,:)=tgradtmp{m};
    fgrad(m,:)=fgradtmp{m};
end;

%% Compare with regular Gabor
g0 = pgauss(L,1);

%[cG,LG] = dgtreal(f,g,a,2*(M-1),L,'timeinv');
[tgradG,fgradG,cG] = gabphasegrad('dgt',f,g0,a,2*(M-1));
tgradG = tgradG(1:M,:); fgradG = fgradG(1:M,:); cG = cG(1:M,:);

b = L/(2*(M-1));

% At this point, tgradG '=' L/2*tgrad and fgradG = fgrad !

reldiff_tgrad = norm(2*tgradG(:)/L-tgrad(:))/norm(tgrad(:));
reldiff_fgrad = norm(fgradG(:)-fgrad(:))/norm(fgrad(:));

% Test correctness of the above assumption by reassigning using tgradG
%tgradGtmp = 2*tgradG/L;
%tgradGtmp = mat2cell(tgradGtmp.',1408,ones(513,1)).';
%
%sr=filterbankreassign(ctmp,tgradGtmp,fgradtmp,a,cfreq);

% Use Gabor version of heapintreal:
PH0=comp_heapintreal(abs(c),tgrad.*L/2,fgrad,a,2*(M-1),1e-4,1);
cG0=abs(c).*exp(1i*PH0);

figure(1);
plotdgtrealphasediff(angle(c),angle(cG0),abs(c),1e-4,a,2*(M-1));

gd0 = gabdual(g0,a,2*(M-1),L);
f_rec0 = idgtreal(cG0,gd0,a,2*(M-1),'timeinv');
soundsc(f_rec0,44100);

%% Uniform Gabor filter bank, heap Integration stuff

PH = comp_heapintreal_ufb(abs(c),tgrad,fgrad,cfreq,a,M,1e-4,1);
c2=abs(c).*exp(1i*PH);

figure(2);
plotdgtrealphasediff(angle(c),angle(c2),abs(c),1e-4,a,2*(M-1));

gd=filterbankrealdual(g,a,L);

f_rec=2*real(ifilterbank(c2.',gd,a));

soundsc(f_rec,44100);

%% Uniform ERB filter bank heap integration stuff

[g,a,cfreq,L] = erbfilters(fs,Ls,'uniform','redmul',2,'spacing',1/6,'bwmul',1/3,'gauss');
cfreq = 2*cfreq/fs;

a = a(1);

[tgradtmp,fgradtmp,~,ctmp] = filterbankphasegrad(f,g,a,L);
 
M = size(ctmp,1);
N = length(ctmp{1});
c=zeros(M,N);
tgrad = c; fgrad = c;
for m=1:M    
    c(m,:)=ctmp{m};    
    tgrad(m,:)=tgradtmp{m};
    fgrad(m,:)=fgradtmp{m};
end;

PH0 = comp_heapintreal_ufb(abs(c),tgrad,fgrad,cfreq,a,M,1e-4,1);
c2=abs(c).*exp(1i*PH0);

figure(3);
plotdgtrealphasediff(angle(c),angle(c2),abs(c),1e-4,a,2*(M-1));

gd=filterbankrealdual(g,a,L);

f_rec=2*real(ifilterbank(c2.',gd,a));

soundsc(f_rec,44100);

%% Same for complex filter bank
[g,a,cfreq,L] = erbfilters(fs,Ls,'uniform','redmul',2,'spacing',1/6,'bwmul',1/3,'complex');
cfreq = 2*cfreq/fs;

a = a(1);

[tgradtmp,fgradtmp,~,ctmp] = filterbankphasegrad(f,g,a,L);
 
M = size(ctmp,1);
N = length(ctmp{1});
c=zeros(M,N);
tgrad = c; fgrad = c;
for m=1:M    
    c(m,:)=ctmp{m};    
    tgrad(m,:)=tgradtmp{m};
    fgrad(m,:)=fgradtmp{m};
end;

PH0 = comp_heapint_ufb(abs(c),tgrad,fgrad,cfreq,a,1e-10,1);
c2=abs(c).*exp(1i*PH0);

figure(4);
plotdgtrealphasediff(angle(c),angle(c2),abs(c),1e-4,a,M);

gd=filterbankdual(g,a,L);

f_rec=real(ifilterbank(c2.',gd,a));

soundsc(f_rec,44100);
