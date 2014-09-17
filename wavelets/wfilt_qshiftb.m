function [h,g,a,info] = wfilt_qshiftb(N)
%WFILT_QSHIFTB  Improved Orthogonality and Symmetry properties 
%
%   Usage: [h,g,a] = wfilt_qshiftb(N);
%
%   `[h,g,a]=wfilt_qshiftb(N)` with $N \in {1,2,3,4,5,6,7}$ returns
%   Kingsbury's Q-shift wavelet filters for tree B.
%
%   Examples:
%   ---------
%   :::
%     figure(1);
%     wfiltinfo('qshiftb3');
% 
%   References: king00 king03
%

% AUTHOR: Zdenek Prusa

[ha,~,a,info] = wfilt_qshifta(N);

hlp = ha{1}.h;
offset = -(numel(hlp)/2 - 1); 
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



