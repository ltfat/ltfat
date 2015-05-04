%DEMO_IMAGECOMPRESSION  Image compression using N-term approximation
%
%   This demo shows how to perform a simple imagecompression using either
%   a Wilson basis or a Wavelet. The compression step is done by
%   retaining only 5% of the coefficients. 
%
%   .. figure::
%
%      Wilson and WMDCT basis
%
%      This right figure shows the image compressed using a |dwilt| basis with
%      8 channels. This corresponds quite closely to JPEG compression,
%      except that the borders between neigbouring blocs are smoother,
%      since the |dwilt| uses a windowing function.
%
%      The left figure shows the same, now
%      using a MDCT basis. The MDCT produces more visible artifacts, as the
%      slowest changing frequency in each block has a half-wave
%      modulation. This is visible on otherwise smooth backgrounds.
%
%   .. figure::
%
%      Wavelet
%
%      The Wavelet used is the DB6 with $J=5$ levels. On the right figure
%      the standard layout has been used, on the left the tensor layout
%      was used.

%% Parameters for the problem formulation

% Use the cameraman image
f=cameraman;

% Ratio to keep
r=0.05;

%% Parameters for the Wilson systems
% Analysis window
ga='itersine';

% synthesis window
gs='itersine';

% No. of channels
M=8;

%% Parameters for the Wavelet system
% Analysis filters
wa='db6';

% Synthesis filters
ws='db6';

% No. of levels 
J=5;

%% Compute the Wilson figures
figure(1);

subplot(1,2,1);
c_dwilt =dwilt2(f,ga,M);
cc_dwilt=largestr(c_dwilt,r);
r_dwilt =idwilt2(cc_dwilt,gs);
imagesc(r_dwilt);
colormap(gray), axis('image');

subplot(1,2,2);
c_wmdct =wmdct2(f,ga,M);
cc_wmdct=largestr(c_wmdct,r);
r_wmdct =iwmdct2(cc_wmdct,gs);
imagesc(r_wmdct);
colormap(gray), axis('image');

%% Compute the Wavelet figures

figure(2);

subplot(1,2,1);
c_fwt =fwt2(f,wa,J);
[cc_fwt,n]=largestr(c_fwt,r);
r_fwt =ifwt2(cc_fwt,ws,J);
imagesc(r_fwt);
colormap(gray), axis('image');

subplot(1,2,2);
c_fwt =fwt2(f,wa,J,'tensor');
cc_fwt=largestr(c_fwt,r);
r_fwt =ifwt2(cc_fwt,ws,J,'tensor');
imagesc(r_fwt);
colormap(gray), axis('image');
