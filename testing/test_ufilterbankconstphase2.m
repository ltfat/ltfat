Ls = 44100;
% f = zeros(Ls,1);
% ind = (1:Ls).'./Ls;
% basefreq = 110;
% fac = 1;
% for kk = 1:7*fac+1
%     f = f + sin(2*pi*2^((kk-1)./fac)*basefreq*ind);
% end
% fs = Ls;
%[f,fs] = gspi;
%f = postpad(f,Ls);
[f,fs] = cocktailparty;
f = f(3*fs + (1:Ls));

%% Uniform ERB filter bank heap integration stuff

spacing = 1/2;
bwmul = 2*spacing;

[g,a,cfreq,L] = erbfilters(fs,Ls,'uniform','redmul',2,'spacing',spacing,'bwmul',bwmul,'gauss');
tfr = getgausstfr(cfreq,fs,L,'erb','bwmul',bwmul);
cfreq = 2*cfreq/fs;

a = a(1);
% 
f = postpad(f,L);
[tgradtmp,fgradtmp,~,ctmp] = filterbankphasegrad(f,g,a,L);
 
M = size(ctmp,1);
N = length(ctmp{1});
c=zeros(N,M);
tgrad = c; fgrad = c;
for m=1:M    
    c(:,m)=ctmp{m};    
    tgrad(:,m)=tgradtmp{m};
    fgrad(:,m)=fgradtmp{m};
end;

[c0,~,~,tgrad0,fgrad0]=ufilterbankconstphase(abs(c),a(1),tfr,cfreq,'real');

% figure(1); close; figure(1); subplot(2,1,1);
% plotfilterbank(abs(tgrad0.'-tgrad.')./abs(tgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg
% 
% subplot(2,1,2);
% plotfilterbank(abs(fgrad0.'-fgrad.')./abs(fgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg

%gd = filterbankrealdual(g,a,L);

figure(2); close; figure(2);
plotdgtrealphasediff(angle(c0),angle(c),abs(c),1e-4,a,M);

%f_rec = 2*real(ifilterbank(c0.',gd,a));
%soundsc(f_rec,44100);

%% Unwarped warped examples

fmax = fs/2;
bins = 1/spacing;

%% ERB alternative
[g,a,cfreq,L]=nonwarpedfilters(@freqtoerb,@erbtofreq,fs,0,fmax,bins,Ls,...
                     'bwmul',bwmul,'real','uniform','gauss');
                 
tfr = getgausstfr_warp(@freqtoerb,@erbtofreq,cfreq,fs,L,'bwmul',bwmul);
cfreq = 2*cfreq/fs;

a = a(1);
% 
f = postpad(f,L);
[tgradtmp,fgradtmp,~,ctmp] = filterbankphasegrad(f,g,a,L);
 
M = size(ctmp,1);
N = length(ctmp{1});
c=zeros(N,M);
tgrad = c; fgrad = c;
for m=1:M    
    c(:,m)=ctmp{m};    
    tgrad(:,m)=tgradtmp{m};
    fgrad(:,m)=fgradtmp{m};
end;

[c0,~,~,tgrad0,fgrad0]=ufilterbankconstphase(abs(c),a(1),tfr,cfreq,'real');

% figure(3); close; figure(3); subplot(2,1,1);
% plotfilterbank(abs(tgrad0.'-tgrad.')./abs(tgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg
% 
% subplot(2,1,2);
% plotfilterbank(abs(fgrad0.'-fgrad.')./abs(fgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg

%gd = filterbankrealdual(g,a,L);

figure(4); close; figure(4);
plotdgtrealphasediff(angle(c0),angle(c),abs(c),1e-4,a,M);

%f_rec = 2*real(ifilterbank(c0.',gd,a));
%soundsc(f_rec,44100);

%% constant-Q
fmin = 50;
warpfun_log = @(x) 10*log(x);
invfun_log = @(x) exp(x/10);
bins = 4
bwmul = 2/bins;


[g,a,cfreq,L]=nonwarpedfilters(warpfun_log,invfun_log,fs,fmin,fmax,bins,Ls,...
                     'bwmul',bwmul,'real','uniform','gauss');
                 
tfr = getgausstfr_warp(warpfun_log,invfun_log,cfreq,fs,L,'bwmul',bwmul);
cfreq = 2*cfreq/fs;

a = a(1);
% 
f = postpad(f,L);
[tgradtmp,fgradtmp,~,ctmp] = filterbankphasegrad(f,g,a,L);
 
M = size(ctmp,1);
N = length(ctmp{1});
c=zeros(N,M);
tgrad = c; fgrad = c;
for m=1:M    
    c(:,m)=ctmp{m};    
    tgrad(:,m)=tgradtmp{m};
    fgrad(:,m)=fgradtmp{m};
end;

[c0,~,~,tgrad0,fgrad0]=ufilterbankconstphase(abs(c),a(1),tfr,cfreq,'real'); 
 
% figure(5); close; figure(5); subplot(2,1,1);
% plotfilterbank(abs(tgrad0.'-tgrad.')./abs(tgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg
% 
% subplot(2,1,2);
% plotfilterbank(abs(fgrad0.'-fgrad.')./abs(fgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg

%gd = filterbankrealdual(g,a,L);

figure(6); close; figure(6);
plotdgtrealphasediff(angle(c0),angle(c),abs(c),1e-4,a,M);

%f_rec = 2*real(ifilterbank(c0.',gd,a));
%soundsc(f_rec,44100);

%% Square root
fmin = 0;
warpfun_sq = @(x) sign(x).*((1+abs(x)/4).^(1/2)-1);
invfun_sq = @(x) 4*sign(x).*((1+abs(x)).^2-1);

[g,a,cfreq,L]=nonwarpedfilters(warpfun_sq,invfun_sq,fs,fmin,fmax,bins,Ls,...
                     'bwmul',bwmul,'real','uniform','gauss');
                 
tfr = getgausstfr_warp(warpfun_sq,invfun_sq,cfreq,fs,L,'bwmul',bwmul);
cfreq = 2*cfreq/fs;

a = a(1);
% 
f = postpad(f,L);
[tgradtmp,fgradtmp,~,ctmp] = filterbankphasegrad(f,g,a,L);
 
M = size(ctmp,1);
N = length(ctmp{1});
c=zeros(N,M);
tgrad = c; fgrad = c;
for m=1:M    
    c(:,m)=ctmp{m};    
    tgrad(:,m)=tgradtmp{m};
    fgrad(:,m)=fgradtmp{m};
end;

[c0,~,~,tgrad0,fgrad0]=ufilterbankconstphase(abs(c),a(1),tfr,cfreq,'real');

% figure(7); close; figure(7); subplot(2,1,1);
% plotfilterbank(abs(tgrad0.'-tgrad.')./abs(tgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg
% 
% subplot(2,1,2);
% plotfilterbank(abs(fgrad0.'-fgrad.')./abs(fgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg

%gd = filterbankrealdual(g,a,L);

figure(8); close; figure(8);
plotdgtrealphasediff(angle(c0),angle(c),abs(c),1e-4,a,M);

%f_rec = 2*real(ifilterbank(c0.',gd,a));
%soundsc(f_rec,44100);

%% Root of order 4                 
warpfun_4 = @(x) 8*sign(x).*((1+abs(x)).^(1/4)-1);
invfun_4 = @(x) sign(x).*((1+abs(x)/8).^4-1);

[g,a,cfreq,L]=nonwarpedfilters(warpfun_4,invfun_4,fs,fmin,fmax,bins,Ls,...
                     'bwmul',bwmul,'real','uniform','gauss');
                 
tfr = getgausstfr_warp(warpfun_4,invfun_4,cfreq,fs,L,'bwmul',bwmul);
cfreq = 2*cfreq/fs;

a = a(1);
% 
f = postpad(f,L);
[tgradtmp,fgradtmp,~,ctmp] = filterbankphasegrad(f,g,a,L);
 
M = size(ctmp,1);
N = length(ctmp{1});
c=zeros(N,M);
tgrad = c; fgrad = c;
for m=1:M    
    c(:,m)=ctmp{m};    
    tgrad(:,m)=tgradtmp{m};
    fgrad(:,m)=fgradtmp{m};
end;

[c0,~,~,tgrad0,fgrad0]=ufilterbankconstphase(abs(c),a(1),tfr,cfreq,'real');

% figure(9); close; figure(9); subplot(2,1,1);
% plotfilterbank(abs(tgrad0.'-tgrad.')./abs(tgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg
% 
% subplot(2,1,2);
% plotfilterbank(abs(fgrad0.'-fgrad.')./abs(fgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,1]); shg

%gd = filterbankrealdual(g,a,L);

figure(10); close; figure(10);
plotdgtrealphasediff(angle(c0),angle(c),abs(c),1e-4,a,M);

%f_rec = 2*real(ifilterbank(c0.',gd,a));
%soundsc(f_rec,44100);
