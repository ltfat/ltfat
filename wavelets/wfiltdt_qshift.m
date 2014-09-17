function [h,g,a,info] = wfiltdt_qshift(N)
%WFILTDT_QSHIFT  Improved Orthogonality and Symmetry properties 
%
%   Usage: [h,g,a] = wfiltdt_qshift(N);
%
%   `[h,g,a]=wfiltdt_qshift(N)` with $N \in {1,2,3,4,5,6,7}$ returns 
%   Kingsbury's Q-shift filters suitable for dual-tree complex wavelet 
%   transform.
%   Filters in both trees are orthogonal and based on a single prototype
%   low-pass filter with a quarter sample delay. Other filters are
%   derived by modulation and time reversal such that they fulfil the
%   half-sample delay difference between the trees.   
%
%   Examples:
%   ---------
%   :::
%     wfiltdtinfo('qshift3');
% 
%   References: king00 king03
%

% AUTHOR: Zdenek Prusa


[h(:,1),g(:,1),a,info] = wfilt_qshifta(N);
[h(:,2),g(:,2)] = wfilt_qshiftb(N);

% Default first and leaf filters
% They are chosen to be orthonormal near-symmetric here in order not to
% break the orthonormality of the overal representation.
[info.defaultfirst, info.defaultfirstinfo] = fwtinit('symorth1');
[info.defaultleaf, info.defaultleafinfo] = ...
    deal(info.defaultfirst,info.defaultfirstinfo);


