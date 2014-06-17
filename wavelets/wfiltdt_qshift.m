function [h,g,a,info] = wfiltdt_qshift(N)
%WFILT_SYMORTH  Improved Orthogonality and Symmetry properties 
%
%   Usage: [h,g,a] = wfilt_symorth(N);
%
%   `[h,g,a]=wfilt_qshift(N)` with *N\in {1,2,3}*
%
%   Examples:
%   ---------
%   :::
%     figure(1);
%     wfiltinfo('ana:symds1');
% 
%   References: kingsbury2000
%

info.istight = 1;
a = [2;2;2;2];

switch(N)
 case 1
    % Example 1. from the reference. Symmetric near-orthogonal
    hlp = [
          0.03516384   % z^4
          0             
         -0.08832942   % z^2    
          0.23389032   % z^1
          0.76027237   % z^0 <-- origin
          0.58751830   % z^-1
          0
         -0.11430184   % z^-3
          0
          0
    ];

    d = 4;

case 2
    % Example 2. From the reference. 
    hlp = [
          0.00325314
         -0.00388321
          0.03466035
         -0.03887280
         -0.11720389
          0.27529538
          0.75614564 % <-- origin
          0.56881042
          0.01186609
         -0.10671180
          0.02382538
          0.01702522
         -0.00543948
         -0.00455690
    ];

    d = 6;
   
case 3
    
      % Example 3. From the reference. 
    hlp = [
        -0.00228413   % z^8
         0.00120989   % z^7  
        -0.01183479   % z^6
         0.00128346   % z^5
         0.04436522   % z^4
        -0.05327611   % z^3
        -0.11330589   % z^2
         0.28090286   % z^1
         0.75281604   % z^0 <-- origin
         0.56580807   % z^-1
         0.02455015   % z^-2
        -0.12018854   % z^-3
         0.01815649   % z^-4
         0.03152638   % z^-5
        -0.00662879   % z^-6
        -0.00257617   % z^-7
         0.00127756   % z^-8
         0.00241187   % z^-9
    ];

    d = 8;
   
  otherwise
        error('%s: No such filters.',upper(mfilename)); 

end
    range = (0:numel(hlp)-1) -d;
    
    % Create the filters according to the reference paper
    harr = [...
            flipud(hlp),...
            (-1).^range.'.*hlp,...
            hlp,...
            flipud((-1).^(range).'.*hlp),...
            ];
        
    % Reverese the filters to obtain the synthesis filters
    harr = flipud(harr);

    % d gets changed
    info.d = [d, d, d, d] +1;

garr = harr;  
h=mat2cell(harr,size(harr,1),ones(1,size(harr,2)));
g=mat2cell(garr,size(garr,1),ones(1,size(garr,2)));


