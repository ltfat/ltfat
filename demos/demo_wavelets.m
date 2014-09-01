%DEMO_WAVELETS  Wavelet filterbanks
%
%   This demo exemplifies the use of the wavelet filterbank trees. All 
%   representations use "least asymmetric" Daubechies wavelet orthonormal
%   filters 'sym8' (8-regular, length 16).
%
%   .. figure::
%
%      DWT representation
%
%      The filterbank tree consists of 11 levels of iterated 2-band basic
%      wavelet filterbank, where only the low-pass output is further 
%      decomposed. This results in 12 bands with octave resolution. 
%
%   .. figure::
%
%      8-band DWT respresentation      
%
%      The filterbank tree (effectively) consists of 3 levels of iterated
%      8-band basic wavelet filterbank resulting in 22 bands. Only the
%      low-pass output is decomposed at each level.
%
%   .. figure::
%
%      Full Wavelet filterbank tree representation
%
%      The filterbank tree depth is 8 and it is fully decomposed meaning
%      both outputs (low-pass and high-pass) of the basic filterbank is
%      plot further. This results in 256 bands linearly covering the 
%      frequency axis. 
%


% Read the test signal and crop it to the range of interest
[f,fs]=gspi;
f=f(10001:100000);
dr=50;
   
    
figure(1);
[c1,info]=fwt(f,'sym8',11);
plotwavelets(c1,info,fs,'dynrange',dr);



figure(2);
[c2,info]=wfbt(f,{'sym8',3,'quadband'});
plotwavelets(c2,info,fs,'dynrange',dr);


figure(3);
[c3,info]=wfbt(f,{'sym8',8,'full'});
plotwavelets(c3,info,fs,'dynrange',dr);


figure(4);
[c4,info]=wfbt(f,{'symorth3',8,'full'});
plotwavelets(c4,info,fs,'dynrange',dr);

figure(5);
[c5,info]=dtwfbreal(f,{'qshift5',8,'full','first','symorth3'});
plotwavelets(c5,info,fs,'dynrange',dr);




