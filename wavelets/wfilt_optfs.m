function [h,g,a,info] = wfilt_optfs(N)
%WFILT_OPTFS Improved frequency separation filters 
%   Usage: [h,g,a] = wfilt_optfs(N);
%
%   `[h,g,a]=wfilt_optfs(N)` returns wavelet filters with better frequency
%   band sepparation (compared to Daubeschies filters of the same length)
%   having at least two vanishing moments.
%
%   References:  paiva2012opt
%
%
a = [2,2];
h = cell(2,1);
switch(N)
    case 1
     h{1}=[0.3452 
           0.7549 
           0.5039 
          -0.0636 
          -0.2117 
           0.0476 
           0.0697 
          -0.0319]; 
    case 2
     h{1}=[0.2657 
           0.7055 
           0.5955 
           0.0292 
          -0.2416
          -0.0248
           0.1238 
          -0.0165 
          -0.0363
           0.0137];
    case 3
     h{1}=[0.2670 
           0.6761 
           0.6096 
           0.0688 
          -0.2564
          -0.0577
           0.1334 
           0.0297 
          -0.0792
           0.0031
           0.0327
          -0.0129];
    case 4
     h{1}=[0.2183 
           0.6298 
           0.6509 
           0.1535 
          -0.2487
          -0.1174
           0.1317 
           0.0687 
          -0.0788
          -0.0345
           0.0534
           0.0002
          -0.0197
           0.0068];
    case 5
     h{1}=[0.2181 
           0.6109 
           0.6541 
           0.1797 
          -0.2430
          -0.1420
           0.1249 
           0.0959 
          -0.0806
          -0.0552
           0.0531
           0.0286
          -0.0378
          -0.0043
           0.0184
          -0.0066];
    case 6
     h{1}=[0.1847 
           0.5703 
           0.6728 
           0.2491 
          -0.2185
          -0.1849
           0.1026 
           0.1253 
          -0.0624
          -0.0823
           0.0455
           0.0486
          -0.0329
          -0.0265
           0.0273
           0.0036
          -0.0120
           0.0039];  
    case 7
     h{1}=[0.1838
           0.5566 
           0.6719 
           0.2682 
          -0.2071
          -0.2012
           0.0905 
           0.1406 
          -0.0529
          -0.0958
           0.0359
           0.0657
          -0.0301
          -0.0400
           0.0240
           0.0220
          -0.0203
          -0.0052
           0.0114
          -0.0038]; 

    otherwise
        error('%s: No such optimized orthonormal wavelet filter.',upper(mfilename));
end;


flen = length(h{1});
h{2} = (-1).^(1:flen).'.*h{1}(end:-1:1);
g = h;
info.istight = 1;
