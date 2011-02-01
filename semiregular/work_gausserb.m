% In this file we construct a uniform filterbank using a Gaussian window for
% each channel. The center frequencies are equidistantly spaced on an ERB-scale,
% and the width of the filter are choosen to match the auditory filter bandwidth as
% determined by Moore.
%
% Each channel is subsampled by a factor of 4 (a=4)
%
% The filterbank cover only the positive frequencies, so we must use filterbankrealdual and
% filterbankrealbounds

f=greasy;
fs=16000;
L=length(f);
a=4;

% Number of channels, slightly less than 1 ERB(Cambridge) per channel.
M=ceil(freqtoerb(fs/2));

% Compute center frequencies.
fc=erbspace(0,fs/2,M);

g=cell(M,1);

% Account for the bandwidth specification to PGAUSS is for 79% of the ERB(integral).  
bw=audfiltbw(fc)/0.79*1.5;

for m=1:M
  g{m}=pgauss(L,'fs',fs,'bw',bw(m),'cf',fc(m));
end;

disp('Frame bound ratio, should be close to 1 if the filters are choosen correctly.');
filterbankrealbounds(g,a)

if 0
gd=filterbankrealdual(g,a);

coef=ufilterbank(f,g,a);
r=2*real(iufilterbank(coef,gd,a));

disp('Error in reconstruction, should be close to zero.');
norm(f-r)
end;
