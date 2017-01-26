Ls = 44100;
[f,fs] = gspi;
f = f(1:Ls);

% Create an UNIFORM filterbank
a = 32;
M = 1024;
tfr = 1/4;

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
plotfilterbank(abs(tgrad0.'-tgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,.01]); shg

subplot(2,1,2);
plotfilterbank(abs(fgrad0.'-fgrad.'),a,'fc',22050*cfreq,'lin','clim',[0,2]); shg



% %% Uniform ERB filter bank heap integration stuff
% 
% [g,a,cfreq,L] = erbfilters(fs,Ls,'uniform','redmul',2,'spacing',1/6,'bwmul',1/3,'gauss');
% cfreq = 2*cfreq/fs;
% 
% a = a(1);
% 
% [tgradtmp,fgradtmp,~,ctmp] = filterbankphasegrad(f,g,a,L);
%  
% M = size(ctmp,1);
% N = length(ctmp{1});
% c=zeros(M,N);
% tgrad = c; fgrad = c;
% for m=1:M    
%     c(m,:)=ctmp{m};    
%     tgrad(m,:)=tgradtmp{m};
%     fgrad(m,:)=fgradtmp{m};
% end;
% 
% PH0 = comp_heapintreal_ufb(abs(c),tgrad,fgrad,cfreq,a,M,1e-4,1);
% c2=abs(c).*exp(1i*PH0);
% 
% figure(3);
% plotdgtrealphasediff(angle(c),angle(c2),abs(c),1e-4,a,2*(M-1));
% 
% gd=filterbankrealdual(g,a,L);
% 
% f_rec=2*real(ifilterbank(c2.',gd,a));
% 
% soundsc(f_rec,44100);
% 
% %% Same for complex filter bank
% [g,a,cfreq,L] = erbfilters(fs,Ls,'uniform','redmul',2,'spacing',1/6,'bwmul',1/3,'complex');
% cfreq = 2*cfreq/fs;
% 
% a = a(1);
% 
% [tgradtmp,fgradtmp,~,ctmp] = filterbankphasegrad(f,g,a,L);
%  
% M = size(ctmp,1);
% N = length(ctmp{1});
% c=zeros(M,N);
% tgrad = c; fgrad = c;
% for m=1:M    
%     c(m,:)=ctmp{m};    
%     tgrad(m,:)=tgradtmp{m};
%     fgrad(m,:)=fgradtmp{m};
% end;
% 
% PH0 = comp_heapint_ufb(abs(c),tgrad,fgrad,cfreq,a,1e-10,1);
% c2=abs(c).*exp(1i*PH0);
% 
% figure(4);
% plotdgtrealphasediff(angle(c),angle(c2),abs(c),1e-4,a,M);
% 
% gd=filterbankdual(g,a,L);
% 
% f_rec=real(ifilterbank(c2.',gd,a));
% 
% soundsc(f_rec,44100);
