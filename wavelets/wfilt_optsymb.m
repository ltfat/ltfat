function [h,g,a,info] = wfilt_optsymb(N)
%WFILT_OPTSYMB  Optimizatized Symmetric Self-Hilbertian Filters 
%
%   Usage: [h,g,a] = wfilt_optsymb(N);
%
%   `[h,g,a]=wfiltdt_optsymb(N)` with $N \in {1,2,3}$ returns filters
%   suitable with optimized symmetry suitable for for dual-tree complex 
%   wavelet transform tree B.
%
%   Examples:
%   ---------
%   :::
%     wfiltinfo('optsymb3');
% 
%   References: dubase08
%

% AUTHOR: Zdenek Prusa


[ha,~,a,info] = wfilt_optsyma(N);


hlp = ha{1}.h;
offset = -(numel(hlp)/2); 
range = (0:numel(hlp)-1) + offset;
    
% Create the filters according to the reference paper.
%
% REMARK: The phase of the alternating +1 and -1 is crucial here.
%         
    harr = [...
            flipud(hlp),...
            (-1).^(range).'.*hlp,...
            ];
        

htmp=mat2cell(harr,size(harr,1),ones(1,size(harr,2)));

h(1:2,1) = cellfun(@(hEl)struct('h',hEl,'offset',offset),htmp(1:2),...
                   'UniformOutput',0);
g = h;
