function [h,g,a,info] = wfiltdt_oddeven(N)
%WFILTDT_ODDEVEN  Kingsbury's symmetric odd and even filters
%
%   Usage: [h,g,a] = wfiltdt_oddeven(N);
%
%   `[h,g,a]=wfilt_oddeven(N)` with *N\in {1}*
%
%   Examples:
%   ---------
%   :::
%     figure(1);
%     wfiltinfo('ana:oddeven1');
% 
%   References: king02
%


a = [2;2];

switch(N)
 case 1
    % Example 1. from the reference. Symmetric near-orthogonal
    garr = [
             0           0           0              0
             0           0           0             -7.062639508928571e-05 
             0          -0.0004645   0              0        
             0           0.0013349  -0.0017578125   1.341901506696429e-03
            -0.0058109   0.0022006   0             -1.883370535714286e-03
             0.0166977  -0.0130127   0.022265625   -7.156808035714285e-03
            -0.0000641   0.0015360  -0.046875       2.385602678571428e-02
            -0.0834914   0.0869008  -0.0482421875   5.564313616071428e-02
             0.0919537   0.0833552   0.2968750     -5.168805803571428e-02
             0.4807151  -0.4885957   0.55546875    -2.997576032366072e-01
             0.4807151   0.4885957   0.2968750      5.594308035714286e-01           
             0.0919537  -0.0833552  -0.0482421875  -2.997576032366072e-01
            -0.0834914  -0.0869008  -0.046875      -5.168805803571428e-02
            -0.0000641  -0.0015360   0.022265625    5.564313616071428e-02
             0.0166977   0.0130127   0              2.385602678571428e-02
            -0.0058109  -0.0022006  -0.0017578125  -7.156808035714285e-03
             0          -0.0013349   0             -1.883370535714286e-03
             0           0.0004645   0              1.341901506696429e-03
             0           0           0              0        
             0           0           0             -7.062639508928571e-05
    ];

    % This scaling is not in the reference paper, but it is here to be
    % consistent
    garr = garr*sqrt(2);
    %garr = normalize(garr,'energy');
    
    offset = -10;

  otherwise
        error('%s: No such filters.',upper(mfilename)); 

end

    %garr = [garr(:,3:4),garr(:,1:2)];
    modrange = (-1).^((0:size(garr,1)-1) + offset).';
    modrange2 = (-1).^((0:size(garr,1)-1) + offset+1).';
    
    harr =       [garr(:,2).*modrange,...
                   garr(:,1).*modrange2,...
                   garr(:,4).*modrange,...
                   garr(:,3).*modrange2
                  ];
            
   
% In the biorthogonal case, the filters do not get time reversed
garr = flipud(garr);
  
htmp=mat2cell(harr,size(harr,1),ones(1,size(harr,2)));
h(1:2,1) = cellfun(@(hEl)struct('h',hEl,'offset',offset),htmp(1:2),...
                   'UniformOutput',0);
h(1:2,2) = cellfun(@(hEl)struct('h',hEl,'offset',offset),htmp(3:4),...
                   'UniformOutput',0);

gtmp=mat2cell(garr,size(garr,1),ones(1,size(garr,2)));

g(1:2,1) = cellfun(@(gEl)struct('h',gEl,'offset',offset),gtmp(1:2),...
                   'UniformOutput',0);
g(1:2,2) = cellfun(@(gEl)struct('h',gEl,'offset',offset),gtmp(3:4),...
                   'UniformOutput',0);

 
[info.defaultfirst, info.defaultfirstinfo] = fwtinit('symorth3');
[info.defaultleaf, info.defaultleafinfo] = ...
    deal(info.defaultfirst,info.defaultfirstinfo);


