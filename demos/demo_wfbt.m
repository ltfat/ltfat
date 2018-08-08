%DEMO_WFBT  Auditory filterbanks built using filterbank tree structures
%
%   This demo shows two specific constructions of wavelet filterbank trees
%   using several M-band filterbanks with possibly different M (making the
%   thee non-homogenous) in order to approximate auditory frequency bands 
%   and a musical scale respectively.  Both wavelet trees produce perfectly
%   reconstructable non-redundant representations. 
%
%   The constructions are the following:
%
%      1) Auditory filterbank, approximately dividing the frequency band 
%         into intervals reminisent of the bask scale.
%
%      2) Musical filterbank, approximately dividing the freq. band into
%         intervals reminicent of the well-tempered musical scale.
%         
%   Shapes of the trees were taken from fig. 8 and fig. 9 from the refernece.
%   Sampling frequency of the test signal is 48kHz as used in the article.
%
%   .. figure:: 
%
%      Frequency responses of the auditory filterbank.
%
%      Both axes are in logarithmic scale.
%   
%   .. figure:: 
%
%      TF plot of the test signal using the auditory filterbank.
% 
%   .. figure:: 
%
%      Frequency responses of the musical filterbank.
%
%      Both axes are in logarithmic scale.
%   
%   .. figure:: 
%
%      TF plot of the test signal using the musical filterbank.
%
%   References: cucla99

% Load a test signal and resample it to 48 kHz since such sampling
% rate is assumed in the reference.
f = dctresample(greasy,3*numel(greasy));
fs = 48000;

% Bark-like filterbank tree
% Creating tree depicted in Figure 8 in the reference. 
w = wfbtinit({'cmband3',1});
w = wfbtput(1,0,'cmband6',w);
w = wfbtput(1,1,'cmband3',w);
w = wfbtput(2,0:1,'cmband5',w);
w = wfbtput(2,2:3,'cmband2',w);

% Convert to filterbank
[g,a] = wfbt2filterbank(w);

% Plot frequency responses
figure(1);
filterbankfreqz(g,a,2*2048,'plot','fs',fs,'dynrange',30,'posfreq','flog');

% Do the transform
[c,info] = wfbt(f,w);
disp('The reconstruction should be close to zero:')
norm(f-iwfbt(c,info))

figure(2);
plotwavelets(c,info,fs,'dynrange',60);

% Well-tempered musical scale filterbank tree
% Creating tree depicted in Figure 9 in the reference. 
w2 = wfbtinit({'cmband4',1});
w2 = wfbtput(1,0:1,'cmband6',w2);
w2 = wfbtput(2,0:1,'cmband4',w2);
w2 = wfbtput(3,1:4,'cmband4',w2);

% Convert to filterbank
[g2,a2] = wfbt2filterbank(w2);
figure(3);
filterbankfreqz(g2,a2,2*2048,'plot','fs',fs,'dynrange',30,'posfreq','flog');


[c2,info2] = wfbt(f,w2);
disp('The reconstruction should be close to zero:')
norm(f-iwfbt(c2,info2))


figure(4);
plotwavelets(c2,info2,fs,'dynrange',60);



