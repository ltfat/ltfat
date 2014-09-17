function [h,g,a,info] = wfiltdt_oddeven(N)
%WFILTDT_ODDEVEN  Kingsbury's symmetric odd and even filters
%
%   Usage: [h,g,a] = wfiltdt_oddeven(N);
%
%   `[h,g,a]=wfilt_oddeven(N)` with $N \in {1}$ returns the original odd
%   and even symmetric filters suitable for dual-tree complex wavelet
%   transform. The filters in individual trees are biorthogonal.
%
%   Examples:
%   ---------
%   :::
%     wfiltdtinfo('ana:oddeven1');
% 
%   References: king02
%

% AUTHOR: Zdenek Prusa


[h(:,1),g(:,1),a,info] = wfilt_oddevena(N);
[h(:,2),g(:,2)] = wfilt_oddevenb(N);
 
[info.defaultfirst, info.defaultfirstinfo] = fwtinit('oddevenb1');
[info.defaultleaf, info.defaultleafinfo] = ...
    deal(info.defaultfirst,info.defaultfirstinfo);


