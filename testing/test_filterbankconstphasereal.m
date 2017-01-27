Ls = 44100;
[f,fs] = gspi;
f = f(1:Ls);

% Create an UNIFORM filterbank
a = 32;
M = 1024;
tfr = 2;

%% Regular uniform Gabor filters

[g,a,M,cfreq,L] = gabrealfilters(a,M,Ls,tfr);
%[g,a,fc] = erbfilters(fs,Ls,'uniform','spacing',1/10,'bwmul',1);

% We will not do no subsampling. 
% This is the main requirement for synchrosqueezing to work.
%a = 1;

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

%% Uniform Gabor filter bank, constructphase stuff
tfr = tfr*ones(size(g))';

[~,~,~,tgrad0,fgrad0]=filterbankconstphasereal(abs(c),g,a(1),tfr,cfreq);

figure(1); close; figure(1); subplot(2,1,1);
plotfilterbank(abs(tgrad0.'-tgrad.')./abs(tgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg

subplot(2,1,2);
plotfilterbank(abs(fgrad0.'-fgrad.')./abs(fgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg

%% Gabor classic

g0 = {'gauss',tfr(1)};

s0 = dgtreal(f,g0,a,2*(M-1),'timeinv');

abss = abs(s0);
logs=log(abss+realmin);
tt=-11;
logs(logs<max(logs(:))+tt)=tt;

difforder = 2;
fgradG = tfr(1)*pderiv(logs,2,difforder)/(2*pi);
% Undo the scaling done by pderiv and scale properly
tgradG = pderiv(logs,1,difforder)/(2*pi*tfr(1))*(2*(M-1)/M);

% Fix the first and last rows .. the
% borders are symmetric so the centered difference is 0
tgradG(1,:) = 0;
tgradG(end,:) = 0;
tgradG = tgradG*2/L;

figure(2); close; figure(2); subplot(2,1,1);
plotfilterbank(abs(tgrad0.'-tgradG.')./abs(tgradG.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg

subplot(2,1,2);
plotfilterbank(abs(fgrad0.'-fgradG.')./abs(fgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg

%% Uniform ERB filter bank heap integration stuff

spacing = 1/6;
bwmul = 2*spacing;

[g,a,cfreq,L] = erbfilters(fs,Ls,'uniform','redmul',2,'spacing',spacing,'bwmul',bwmul,'gauss');
tfr = getgausstfr(cfreq,fs,L,'erb','bwmul',bwmul);
cfreq = 2*cfreq/fs;

a = a(1);

f = postpad(f,L);
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

[c0,~,~,tgrad0,fgrad0]=filterbankconstphasereal(abs(c),g,a(1),tfr,cfreq);

figure(3); close; figure(3); subplot(2,1,1);
plotfilterbank(abs(tgrad0.'-tgrad.')./abs(tgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg

subplot(2,1,2);
plotfilterbank(abs(fgrad0.'-fgrad.')./abs(fgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg

gd = filterbankrealdual(g,a,L);

figure(4); close; figure(4);
plotdgtrealphasediff(angle(c0),angle(c),abs(c),1e-4,a,M);

f_rec = 2*real(ifilterbank(c0.',gd,a));
soundsc(f_rec,44100);

%% Same for complex filter bank
[g,a,cfreq,L] = erbfilters(fs,Ls,'uniform','redmul',2,'spacing',spacing,'bwmul',bwmul,'complex','gauss');
tfr = getgausstfr(abs(cfreq),fs,L,'erb','bwmul',bwmul);
cfreq = 2*cfreq/fs;

a = a(1);

f = postpad(f,L);
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

[c0,~,~,tgrad0,fgrad0]=filterbankconstphase(abs(c),g,a(1),tfr,cfreq);

figure(5); close; figure(5); subplot(2,1,1);
plotfilterbank(abs(tgrad0.'-tgrad.')./abs(tgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg

subplot(2,1,2);
plotfilterbank((fgrad0.'-fgrad.')./abs(fgrad.'),a,'fc',22050*cfreq,'lin','clim',[-1,1]); shg

gd = filterbankdual(g,a,L);

figure(6); close; figure(6);
plotdgtrealphasediff(angle(c0),angle(c),abs(c),1e-4,a,2*(M-1));

f_rec = real(ifilterbank(c0.',gd,a));
soundsc(f_rec,44100);
