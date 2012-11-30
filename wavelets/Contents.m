% LTFAT - Wavelets
%
%   Basic wavelet analysis/synthesis interface functions
%      FWT
%      IFWT
%      WAVELETFB
%      WTFFT
%      IWTFFT
%      
%
%   Plots
%      PLOTWAVC          - Plot wavelet coefficients
%      FREQZFB           - Calculate or plot filterbank frequency responses
%  
%
%
%   Auxilary
%      MULTID            - Creates equivalent one-level multirate identity filterbank
%      WAVFUN            - Aproximate of the continuous scaling and wavelet functions
%      CELL2PACK         - Changes wavelet coefficient storing format
%      PACK2CELL         - Changes wavelet coefficient storing format
%
%   Wavelet filters defined in time-domain
%      WFILT_DB          - DauBechies orthogonal filters (ortonormal base)
%      WFILT_DDEN        - Double-DENsity dwt filters (tight frame)
%      WFILT_DGRID       - Dense GRID framelets (tight frame, symmetric)
%      WFILT_DTREE       - Dual-TREE complex wavelet transform filters (two orthonormal bases)
%      WFILT_HDEN        - Higher DENsity dwt filters (tight frame, frame)  
%      WFILT_OPTFS       - Optimized orthogonal filters with improved Frequency Selectivity (ortonormal base)
%      WFILT_SYMDS       - SYMmetric wavelet Dyadic Siblings (frames)
%
%   Wavelet filters defined in frequency-domain
%      WFREQ_RADWT       - Overcomplete Rational-Dilation wavelet transform
%
%   Continuous wavelets (no perfect reconstruction)