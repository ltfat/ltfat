%DEMO_FILTERBANKS  CQT and ERBLET filterbanks
%
%   This demo shows CQT (Constant Quality Transform) and ERBLET (Equivalent
%   Rectangular Bandwidth -let transform) representations acting as
%   filterbanks  with high and low redundancies. Filterbanks are build such
%   that the painless condition is always satified. Real input signal and 
%   filters covering only the positive frequency range are used. The 
%   redundancy is calculated as a ratio of the number of (complex) 
%   coefficients and the input length times two to account for the storage
%   requirements of complex numbers.
%
%     * The high redundancy representation uses 'uniform' subsampling i.e.
%       all channels are subsampled with the same subsampling factor which
%       is the lowest from the filters according to the painless condition
%       rounded towards zero.
%
%     * The low redundancy representation uses 'fractional' subsampling
%       which results in the least redundant representation still
%       satisfying the painless condition. Actual time positions of atoms 
%       can be non-integer, hence the word fractional.
%
%   .. figure::
%
%      ERBLET representations
%
%      The high-redundany plot (top) consists of 400 channels (~9 filters 
%      per ERB) and low-redundany plot (bottom) consists of 44 channels 
%      (1 filter per ERB).
%
%   .. figure::
%
%      CQT representations
%
%      Both representations consist of 280 channels (32 channels per octave,
%      frequency range 50Hz-20kHz). The high-redundany represention is on 
%      the top and the low-redundancy repr. is on the bottom.
%
%   See also: erbfilters, cqtfilters

% Read the test signal and crop it to the range of interest
[f,fs]=gspi;
f=f(10001:100000);
dr=50;

figure(1);
subplot(2,1,1);

% Create the ERB filterbank using 400 filters linearly spaced in the 
% ERB scale and using uniform subsampling.
[g,a,fc]=erbfilters(fs,numel(f),'M',400,'uniform');

% Compute the filterbank response
c1=filterbank(f,g,a);

% Compute redundancy
erb1_redundancy = 2*sum(1./a);

% Plot the representation
plotfilterbank(c1,a,fc,fs,dr,'audtick');

subplot(2,1,2);

% Create the ERB filterbank using 44 filters linearly spaced in the 
% ERB scale and using fractional subsampling.
[g,a,fc]=erbfilters(fs,numel(f),'fractional');

% Compute the filterbank response
c2=filterbank(f,g,a);

% Compute redundancy
erb2_redundancy = 2*sum(a(:,2)./a(:,1));

% Plot the representation
plotfilterbank(c2,a,fc,fs,dr,'audtick');

fprintf('ERBLET high redundancy %.2f, low redundany %.2f.\n',...
         erb1_redundancy,erb2_redundancy);

figure(2);
subplot(2,1,1);

% Create the CQT filterbank using 32 channels per octave in frequency 
% range 50Hz-20kHz using uniform subsampling.
[g,a,fc] = cqtfilters(fs,50,20000,32,numel(f),'uniform');

% Compute the filterbank response
c3=filterbank(f,g,a);

% Compute redundancy
cqt1_redundancy = 2*sum(1./a);

% Plot the representation
plotfilterbank(c3,a,fc,fs,'dynrange',dr);

subplot(2,1,2);

% Create the CQT filterbank using 32 channels per octave in frequency 
% range 50Hz-20kHz using fractional subsampling.
[g,a,fc] = cqtfilters(fs,50,20000,32,numel(f),'fractional');

% Compute the filterbank response
c4=filterbank(f,g,a);

% Compute redundancy
cqt2_redundancy = 2*sum(a(:,2)./a(:,1));

% Plot the representation
plotfilterbank(c4,a,fc,fs,'dynrange',dr);

fprintf('CQT high redundancy %.2f, low redundany %.2f.\n',...
        cqt1_redundancy,cqt2_redundancy);

