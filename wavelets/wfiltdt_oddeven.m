function [h,g,a,info] = wfiltdt_oddeven(N)
%WFILT_SYMORTH  Kingsbury's symmetric (nearly orthogonal) filters
%
%   Usage: [h,g,a] = wfilt_symorth(N);
%
%   `[h,g,a]=wfilt_symorth(N)` with *N\in {1,2,3}*
%
%   Examples:
%   ---------
%   :::
%     figure(1);
%     wfiltinfo('ana:symds1');
% 
%   References: kingsbury2000
%


a = [2;2;2;2];

switch(N)
 case 1
    % Example 1. from the reference. Symmetric near-orthogonal
    garr = [
           0          -0.0000706   0            0   
           0           0           0           -0.0004645 
          -0.0017581   0.0013419   0            0.0013349  
           0          -0.0018834  -0.00581090   0.0022006 
           0.0222656  -0.0071568   0.01669770  -0.0130127 
          -0.0468750   0.0238560  -0.00006410   0.0015360 
          -0.0482422   0.0556431  -0.08349140   0.0869008 
           0.2968750  -0.0516881   0.09195370   0.0833552 
           0.5554688  -0.2997576   0.48071510  -0.4885957 
           0.2968750   0.5594308   0.48071510   0.4885957             
          -0.0482422  -0.2997576   0.09195370  -0.0833552
          -0.0468750  -0.0516881  -0.08349140  -0.0869008
           0.0222656   0.0556431  -0.00006410  -0.0015360
           0           0.0238560   0.01669770   0.0130127
          -0.0017581  -0.0071568  -0.00581090  -0.0022006
           0          -0.0018834   0           -0.0013349
           0           0.0013419   0            0.0004645 
           0           0           0            0
           0          -0.0000706   0            0 
    ];

    garr = normalize(garr,'energy');

    d = 8;

  otherwise
        error('%s: No such filters.',upper(mfilename)); 

end
    range = (1:size(garr,1)) - d;
    garr = [garr(:,3:4),garr(:,[1,2])];
    
    harr = bsxfun(@times,(-1).^range.',...
                  [garr(:,2),...
                   garr(:,1),...
                   garr(:,4),...
                   garr(:,3)
                  ]);
            
    harr = flipud(harr);
    garr = flipud(garr);
     
info.d = [d,d,d,d];

h=mat2cell(harr,size(harr,1),ones(1,size(harr,2)));
g=mat2cell(garr,size(garr,1),ones(1,size(garr,2)));

info.defaultfirst = 'syn:symorth3'; 


