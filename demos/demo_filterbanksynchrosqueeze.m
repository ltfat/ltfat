%DEMO_FILTERBANKSYNCHROSQUEEZE Filterbank synchrosqueezing and inversion
%
%   The demo shows that the synchrosqueezed filterbank representation can be 
%   directly used to reconstruct the original signal.
%   Since we do not work with a filterbank which forms a tight frame 
%   (its |filterbankresponse| is not constant) the direct reconstruction 
%   (mere summing all the channels) does not work well. We can fix that by
%   filtering (equalizing) the result by the inverse of the overall analysis 
%   filterbank frequency response.
%
%   .. figure::
%      
%      ERBlet spectrogram (top) and synchrosqueezed ERBlet spectrogram (bottom)
%
%      The signal used is the first second from |gspi|. Only the energy of
%      the coefficients is show. Both representations are in fact complex and
%      invertible.
%
%   .. figure::
%
%      Errors of the direct and the equalized reconstructions
%       
%      There is still a small DC offset of the signal obtained by the direct
%      summation.
%
%   References: ltfatnote041

% Get one second of gspi
Ls = 44100;
[f,fs] = gspi;
f = f(1:Ls);

% Create an UNIFORM filterbank
[g,a,fc] = erbfilters(fs,Ls,'uniform','M',400);

% We will not do no subsampling. 
% This is the main requirement for synchrosqueezing to work.
a = 1;

% Compute the time phase derivative and the coefficients
[tgrad,~,~,c] = filterbankphasegrad(f,g,a);

% Do the synchrosqueezing
cs = filterbanksynchrosqueeze(c,tgrad,cent_freqs(fs,fc));

% Plot spectrograms
figure(1);clf;
subplot(2,1,1);
plotfilterbank(c,a,fc,fs,60);
subplot(2,1,2);
plotfilterbank(cs,a,fc,fs,60);

% Reformat the coefficients to matrices
cmat = cell2mat(c.');
csmat = cell2mat(cs.');

% Compute the overall analysis filterbank response.
Frespall = sum(filterbankfreqz(g,a,Ls),2);

% The fiterbank is defined only for positive frequencies, we 
% sum the response with its involution
Frespfull = Frespall + involute(Frespall);

% The filterbank is not tight, so the direct reconstruction 
% will not give a good reconstruction.
% We do that anyway for comparison.
% 
% Just scale so that the response is around 1
C = mean(abs(Frespfull));
Frespfull = Frespfull/C;

% Direct reconstruction from the original coefficients 
fhat1  = 2*real(sum(cmat,2))/C;
% Direct reconstruction from the synchrosqueezed coefficents
fhat1s = 2*real(sum(csmat,2))/C;

% Reconstruction errors
err1 = norm(f-fhat1)/norm(f);
err1s = norm(f-fhat1s)/norm(f);

% We "equalize" the reconstruction by inverse of Frespfull
fhat2 = real(ifft(fft(fhat1)./Frespfull));
fhat2s = real(ifft(fft(fhat1s)./Frespfull));

% Compute errors
err2 = norm(f-fhat2)/norm(f);
err2s = norm(f-fhat2s)/norm(f);

% Plot errors for comparison
figure(2);clf;
title('Error of the reconstruction');
plot([f-fhat1s, f-fhat2s]);
legend({'notequalized','equalized'});

fprintf(['Direct reconstruction MSE:\n   From the coefficients: %e\n',...
         '   From the synchrosqueezed coefficients: %e\n'],err1,err1s);

fprintf(['Equalized reconstruction MSE:\n   From the coefficients: %e\n',...
         '   From the synchrosqueezed coefficients: %e\n'],err2,err2s);

