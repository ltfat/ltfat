function [h,g,a,info] = wfilt_algmband(K)
%WFILT_ALGMBAND  An ALGebraic construction of orthonormal M-BAND wavelets with perfect reconstruction
%   Usage: [h,g,a] = wfilt_algmband(K);
%
%   `[h,g,a]=wfilt_algmband(K)` with $K \in {1,2}$ returns wavelet filters
%   from the reference paper. The filters are 3-band ($K==1$) and 4-band 
%   ($K==2$) with critical subsampling.
%
%   Examples:
%   ---------
%   :::
%     wfiltinfo('algmband1');  
%
%   :::
%     wfiltinfo('algmband2');   
%
%   References:  lin2006algebraic  

% AUTHOR: Zdenek Prusa


switch(K)
   case 1
   % from the paper Example 1.
      garr = [
              0.33838609728386 -0.11737701613483 0.40363686892892
              0.53083618701374 0.54433105395181 -0.62853936105471
              0.72328627674361 -0.01870574735313 0.46060475252131
              0.23896417190576 -0.69911956479289 -0.40363686892892
              0.04651408217589 -0.13608276348796 -0.07856742013185
             -0.14593600755399 0.42695403781698 0.24650202866523
             ];
      a= [3;3;3];
      offset = [-3,-3,-3];
   case 2
      garr = [
              0.0857130200  -0.1045086525 0.2560950163  0.1839986022
              0.1931394393  0.1183282069  -0.2048089157 -0.6622893130
              0.3491805097  -0.1011065044 -0.2503433230 0.6880085746
              0.5616494215  -0.0115563891 -0.2484277272 -0.1379502447
              0.4955029828  0.6005913823  0.4477496752  0.0446493766
              0.4145647737  -0.2550401616 0.0010274000  -0.0823301969
              0.2190308939  -0.4264277361 -0.0621881917 -0.0923899104
             -0.1145361261 -0.0827398180 0.5562313118  -0.0233349758
             -0.0952930728 0.0722022649  -0.2245618041 0.0290655661
             -0.1306948909 0.2684936992  -0.3300536827 0.0702950474
             -0.0827496793 0.1691549718  -0.2088643503 0.0443561794
              0.0719795354  -0.4437039320 0.2202951830  -0.0918374833
              0.0140770701  0.0849964877  0.0207171125  0.0128845052
              0.0229906779  0.1388163056  0.0338351983  0.0210429802
              0.0145382757  0.0877812188  0.0213958651  0.0133066389
             -0.0190928308 -0.1152813433 -0.0280987676 -0.0174753464
             ];
       a= [4;4;4;4];
       offset = [-12,-8,-8,-12];
  otherwise
        error('%s: No such orthonormal M-band wavelet filter bank.',upper(mfilename));
end

g=mat2cell(flipud(garr),size(garr,1),ones(1,size(garr,2)));
g = cellfun(@(gEl,offEl) struct('h',gEl,'offset',offEl),g,num2cell(offset),...
            'UniformOutput',0);

h = g;
info.istight=1;
