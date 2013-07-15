function [h,g,a,info] = wfilt_dgrid(N)
%WFILT_DGRID  Dense GRID framelets (tight frame, symmetric)
%   Usage: [h,g,a] = wfilt_dgrid(N);
%
%   `[h,g,a]=wfilt_dgrid(N)` computes Dense GRID framelets. Redundancy
%   equal to 2.
%
%   References: abdelnour2007dense
%
%



switch(N)
case 1
harr = [
  0              0               0              0   
  -sqrt(2)/32    -sqrt(2)/32     0.0104677975   0
  0              -2*sqrt(2)/32   0.0370823430   0.0104677975
  9*sqrt(2)/32   7*sqrt(2)/32    0.0651503417   0.0370823430
  16*sqrt(2)/32  0               0.0897493772   0.0651503417
  9*sqrt(2)/32   -7*sqrt(2)/32   -0.575618139   0.0897493772
  0              2*sqrt(2)/32    0.3731682798   -0.575618139
  -sqrt(2)/32    sqrt(2)/32      0              0.3731682798
];
case 2
harr = [
  -5*sqrt(2)/256       -5*sqrt(2)/256     0.0422028267     0
  -7*sqrt(2)/256       23*sqrt(2)/256     0.0784808462     0.0422028267 
  35*sqrt(2)/256       -13*sqrt(2)/256    0.0274495253     0.0784808462
  105*sqrt(2)/256      -41*sqrt(2)/256    -0.1272844093    0.0274495253
  105*sqrt(2)/256      41*sqrt(2)/256     -0.4611848140    -0.1272844093   
  35*sqrt(2)/256       13*sqrt(2)/256     0.5488035631     -0.4611848140            
  -7*sqrt(2)/256       -23*sqrt(2)/256    -0.1084675382    0.5488035631 
  -5*sqrt(2)/256       5*sqrt(2)/256      0                -0.1084675382
]; 
case 3
harr = [
  -7*sqrt(2)/1024       -7*sqrt(2)/1024     0.0019452732     0
  -27*sqrt(2)/1024       43*sqrt(2)/1024    -0.0020062621    0.0019452732  
  0                    -80*sqrt(2)/1024     0.0070362139     -0.0020062621
  168*sqrt(2)/1024      8*sqrt(2)/1024      -0.0305577537    0.0070362139
  378*sqrt(2)/1024      138*sqrt(2)/1024    -0.1131305218    -0.0305577537  
  378*sqrt(2)/1024      -138*sqrt(2)/1024   -0.0700905154    -0.1131305218  
  168*sqrt(2)/1024      -8*sqrt(2)/1024     0.0845961181     -0.0700905154  
  0                     80*sqrt(2)/1024     0.6026545312     0.0845961181
  -27*sqrt(2)/1024      -43*sqrt(2)/1024    -0.4804470835    0.6026545312
  -7*sqrt(2)/1024       7*sqrt(2)/1024      0                -0.4804470835
]; 

    otherwise
        error('%s: No such Dense Grid Framelet filters.',upper(mfilename));
end;

garr = flipud(harr);
g=mat2cell(garr,size(garr,1),ones(1,size(garr,2)));
h = g;
a = [2;2;2;2];
info.istight=1;