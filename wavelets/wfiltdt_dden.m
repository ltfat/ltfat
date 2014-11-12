function [h,g,a,info] = wfiltdt_dden(N)
%WFILTDT_DDEN  Double-Density Dual-Tree DWT filters 
%
%   Usage: [h,g,a] = wfiltdt_dden(N);
%
%   `[h,g,a]=wfiltdt_dden(N)` with $N \in {1,2}$ returns filters suitable
%   for dual-tree double density complex wavelet transform. 
%
%   Examples:
%   ---------
%   :::
%     wfiltdtinfo('dden1');
%
%   :::
%     wfiltdtinfo('dden2');
% 
%   References: se04
%

% AUTHOR: Zdenek Prusa

[h(:,1),g(:,1),a,info] = wfilt_ddena(N);
[h(:,2),g(:,2)] = wfilt_ddenb(N);

[info.defaultfirst, info.defaultfirstinfo] = fwtinit('symdden2');
[info.defaultleaf, info.defaultleafinfo] = ...
    deal(info.defaultfirst,info.defaultfirstinfo);




