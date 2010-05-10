%DEMO_PARTIALS  Find partials in speech signal
%
%   This script shows a simple method for finding partials in a speech
%   signal. Partials are peaks in a local spectrum or spectogram, and
%   can be grouped into formants. The latter are linked to resonant
%   frequencies of the related physical system. They are most commonly
%   used in phonetics, searching for the formants of the human vocal tract.
%   This is for demo used for speaker recognition.
%
%R   carmonmultiridge1 fant1 stevensphon1

% AUTHORs:  Peter Balazs
%           Peter Soendergaard 

disp('Type "help demo_partials" to see a description of how this demo works.');

a=64; % hop size
M=1024; % number of frequency bins
freq_resol = 4; % increase frequency resolution of the Gabor windows

SR = 16000; % sampling rate

% Define the input
demoname='greasy';
f=greasy;
fs=length(f);

% Length of transform to do
L=6144;

% Window function
winname='Gauss';
g=pgauss(L,freq_resol*a*M/L);

% Compute the cofficients.
c=dgt(f,g,a,M);

% Get the number of time shifts.
N=size(c,2);

% For voval tract formants four are often used for speaker recognition.
F=5; 

% This should be vectorized.
partials = zeros(F,N);
for n=1:N;
  [partials(:,n),Y] = findmaxima(abs(c(1:M/2,n)),F);
end;

figure(1);
specto = 20*log10(abs(c(1:M/2,:))); % log plot, more related to perception
maxlV = max(max(specto));
minlV = maxlV - 50; % 50 db range 
imagesc(specto,[minlV,maxlV]);
title(sprintf('Log. Amplitude STFT of %s \n Parameters: %s window, Nwin = %g, a = %g, M = %g',demoname,winname,N,a,M));
set(gca,'Ydir','normal');

ylabel('kHz');
yt = get(gca,'ytick');
yt = yt*(SR/M)/1000; % rescaling of y axis for physical units 
 
ytext = sprintf('%-2.1f',yt(1));
for ii=2:length(yt)
     ytext = sprintf('%s | %-2.1f',ytext,yt(ii));
end
set(gca,'yticklabel',ytext);
xlabel('s');
xt = get(gca,'xtick');
xt = xt*a/SR; % rescaling of x axis for physical units
xtext = sprintf('%-2.1f',xt(1));
for ii=2:length(xt)
     xtext = sprintf('%s | %-2.2f',xtext,xt(ii));
end
set(gca,'xticklabel',xtext);

hold on;
partials = partials.';
plot(partials(:,1),'k.:','Linewidth',2);
plot(partials(:,2),'g.:','Linewidth',2);
plot(partials(:,3),'r.:','Linewidth',2);
plot(partials(:,4),'b.:','Linewidth',2);
plot(partials(:,5),'m.:','Linewidth',2);
legend('partial1','partial2','partial3','partial4','partial5');
hold off;

% Formant visualization: Changing the time-frequency resolution, hop size
% and number of frequency bins

a=8; % hop size
M=256; % number of frequency bins
freq_resol = 1/4; % lower frequency resolution of the Gabor windows

% Window function
winname='Gauss';
g=pgauss(L,freq_resol*a*M/L);

% Compute the cofficients.
c=dgt(f,g,a,M);

% Get the number of time shifts.
N=size(c,2);

% display
figure(2);
specto = 20*log10(abs(c(1:M/2,:))); % log plot, more related to perception
maxlV = max(max(specto));
minlV = maxlV - 50; % 50 db range 
imagesc(specto,[minlV,maxlV]);
title(sprintf('Log. Amplitude STFT of %s \n Parameters: %s window, Nwin = %g, a = %g, M = %g',demoname,winname,N,a,M));
set(gca,'Ydir','normal');

ylabel('kHz');
yt = get(gca,'ytick');
yt = yt*(SR/M)/1000; % rescaling of y axis for physical units 
 
ytext = sprintf('%-2.1f',yt(1));
for ii=2:length(yt)
     ytext = sprintf('%s | %-2.1f',ytext,yt(ii));
end
set(gca,'yticklabel',ytext);
xlabel('s');
xt = get(gca,'xtick');
xt = xt*a/SR; % rescaling of x axis for physical units
xtext = sprintf('%-2.1f',xt(1));
for ii=2:length(xt)
     xtext = sprintf('%s | %-2.2f',xtext,xt(ii));
end
set(gca,'xticklabel',xtext);
