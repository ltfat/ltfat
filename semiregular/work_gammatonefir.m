% In this file we construct a uniform filterbank using a the impulse response of a
% 4th order gammatone  for each channel. The center frequencies are equidistantly spaced 
% on an ERB-scale,
% and the width of the filter are choosen to match the auditory filter bandwidth as
% determined by Moore.
%
% Each channel is subsampled by a factor of 4 (a=4)
%
% The filterbank cover only the positive frequencies, so we must use filterbankrealdual and
% filterbankrealbounds

if 0
  % Use greasy, low sampling rate, short signal
  f=greasy;
  fs=16000;
  a=16;
  erbs_per_channel=1;
  filterlength=ceil(length(f)/a)*a;
else
  % glockenspiel, high sampling rate, longer signal
  f=gspi;
  fs=44100;
  a=16;
  erbs_per_channel=1;
  filterlength=5000;
end;

% Determine minimal transform length
Ls=length(f);
L=ceil(filterlength/a)*a;

% Number of channels, slightly less than 1 ERB(Cambridge) per channel.
M=ceil(freqtoerb(fs/2)/erbs_per_channel);

% Compute center frequencies.
fc=erbspace(0,fs/2,M);

g=gammatonefir(fc,fs,filterlength);

% In production code, it is not necessary to call 'filterbankrealbounds',
% this is just for veryfying the setup.
disp('Frame bound ratio, should be close to 1 if the filters are choosen correctly.');
filterbankrealbounds(g,a,L)

% Create reconstruction filters
gd=filterbankrealdual(g,a,L);

coef=ufilterbank(f,g,a);
r=2*real(iufilterbank(coef,gd,a,Ls));

disp('Error in reconstruction, should be close to zero.');
norm(f-r)
  
