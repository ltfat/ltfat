function [h,g,a,info] = wfilt_matlabwtwrapper(wname)
%WFILT_MATLABWTWRAPPER Wrapper of the Matlab Wavelet Toolbox wfilters function
%   Usage: [h,g,a] = wfilt_matlabwtwrapper(wname);
%
%   `[h,g,a]=wfilt_matlabwtwrapper(wname)` calls Matlab Wavelet Toolbox
%   function `wfilters` and passes the parameter *wname*. This function
%   requires the Matlab Wavelet Toolbox.
%

a = [2;2];
[lo,hi,lo_s,hi_s] = wfilters(wname);

h=cell(2,1);
h{1} = lo(:);
h{2} = hi(:);

g=cell(2,1);
g{1} = flipud(lo_s(:));
g{2} = flipud(hi_s(:));

if all(h{1}==g{1}) && all(h{2}==g{2})
  info.istight = 1;
else
  info.istight = 0; 
end

