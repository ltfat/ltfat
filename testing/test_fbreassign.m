function test_fbreassign

L = 44100;

%f = exp(2*pi*1i*((0:44099)/168))';
f = sin(2*pi*((0:44099)/35+((0:44099)/300).^2)) + ...
    sin(2*pi*((0:44099)/10+((0:44099)/300).^2)) + ...
    sin(2*pi*((0:44099)/5-((0:44099)/450).^2));
f = 0.7*f';

%f = zeros(L,1);
%f(10000) = 1;

% wavwrite(f,44100,16,'testclean.wav');
% f = gspi;
% f = f(44101:88200);
% fn = rand(44100,1)-.5;
% fn = fn/norm(fn)*norm(f);
% f = f+fn;
% wavwrite(f,44100,16,'testnoisy.wav');
% f = gspi;
% f = f(44101:88200);

fs = 44100;

[g,a,fc,L]=erbfilters(fs,44100,'fractional','spacing',1/12,'warped','complex');
if length(a) < length(g)
    a = [a;a(end-1:-1:2,:)];
end
tic; [tgrad,fgrad,cs0,c]=filterbankphasegrad(f,g,a,L); PGtime = toc

tic; fc = cent_freqs(fs,fc); CFtime = toc
tic; sr0=filterbankreassign(cs0,tgrad,fgrad,a,fc); RAtime = toc
figure(1); clf;
subplot(211);
plotfilterbank(cs0,a,'fc',fs/2*fc,'db','dynrange',60);
title('ERBlet spectrogram of 3 chirps');
subplot(212);  plotfilterbank(sr0,a,'fc',fs/2*fc,'db','dynrange',60);
title('Reassigned ERBlet spectrogram of 3 chirps');
colormap(flipud(gray));



%% 
clear g g2 g3
Lg = 28;
a0 = 7*ones(168,1);
cfreq0 = [0:84,-83:-1].'/84;
gg = fftshift(firwin('hann',Lg));

for kk = 0:167
g{kk+1}.h = gg;
g{kk+1}.fc = modcent(kk/84,2);
g{kk+1}.offset = -Lg/2; 
g{kk+1}.realonly = 0; 
end

tic; [tgrad,fgrad,c_s]=filterbankphasegrad(f,g,a0,L); PGtimeFIR = toc
tic; sr=filterbankreassign(c_s,tgrad,fgrad,a0,cfreq0); RAtimeFIR = toc

Lg = 882;
a = 90*ones(882,1);
cfreq1 = [0:441,-440:-1].'/441;
gg = fftshift(firwin('hann',Lg));%.*exp(-2*pi*1i*100*(-Lg/2:Lg/2-1).'./L);
for kk = 0:881
g2{kk+1}.H = gg;
g2{kk+1}.L = L;
g2{kk+1}.foff = kk*L/882-Lg/2;
g2{kk+1}.realonly = 0; 
g3{kk+1}.L = L;
g3{kk+1}.H = comp_transferfunction(g2{kk+1},L);
g3{kk+1}.foff = 0;
g3{kk+1}.realonly = 0; 
end

figure(2); clf; subplot(311); 
plotfilterbank(sr,a0,'fc',fs/2*cfreq0,'linabs');

tic; [tgrad,fgrad,c_s2]=filterbankphasegrad(f,g2,a,L); PGtimeBL = toc
tic; sr2=filterbankreassign(c_s2,tgrad,fgrad,a,cfreq1); RAtimeBL = toc
figure(2); subplot(312); plotfilterbank(sr2,a,'fc',fs/2*cfreq1,'linabs');

tic; [tgrad,fgrad,c_s3]=filterbankphasegrad(f,g3,a,L); PGtimeL = toc
tic; sr3=filterbankreassign(c_s3,tgrad,fgrad,a,cfreq1); RAtimeL = toc
figure(2); subplot(313); plotfilterbank(sr3,a,'fc',fs/2*cfreq1,'linabs');
% 
% clear all