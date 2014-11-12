function [h,g,a,info] = wfiltdt_optsym(N)
%WFILTDT_OPTSYM  Optimizatized Symmetric Self-Hilbertian Filters 
%
%   Usage: [h,g,a] = wfiltdt_optsym(N);
%
%   `[h,g,a]=wfiltdt_optsym(N)` with $N \in {1,2,3}$ returns filters
%   suitable for dual-tree complex wavelet transform with optimized 
%   symmetry.
%
%   Examples:
%   ---------
%   :::
%     wfiltdtinfo('optsym3');
% 
%   References: dubase08
%

% AUTHOR: Zdenek Prusa

               
[h(:,1),g(:,1),a,info] = wfilt_optsyma(N);
[h(:,2),g(:,2)] = wfilt_optsymb(N);

% Default first and leaf filters
% They are chosen to be orthonormal near-symmetric here in order not to
% break the orthonormality of the overal representation.
[info.defaultfirst, info.defaultfirstinfo] = fwtinit('symorth1');
[info.defaultleaf, info.defaultleafinfo] = ...
    deal(info.defaultfirst,info.defaultfirstinfo);


