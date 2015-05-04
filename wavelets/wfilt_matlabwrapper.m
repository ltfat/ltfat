function [h,g,a,info] = wfilt_matlabwrapper(wname)
%WFILT_MATLABWRAPPER Wrapper of the Matlab Wavelet Toolbox wfilters function
%   Usage: [h,g,a] = wfilt_matlabwrapper(wname);
%
%   `[h,g,a]=wfilt_matlabwrapper(wname)` calls Matlab Wavelet Toolbox
%   function `wfilters` and passes the parameter *wname* to it. 
%
%   This function requires the Matlab Wavelet Toolbox.
%

% AUTHOR: Zdenek Prusa

if ~exist('wfilters',2)
    error('%s: Matlab Wavelet Toolbox is not present.',upper(mfilename));
end

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

g = cellfun(@(gEl) struct('h',gEl(:),'offset',-numel(gEl)/2),g,'UniformOutput',0);
h = cellfun(@(hEl) struct('h',hEl(:),'offset',-numel(hEl)/2),h,'UniformOutput',0);
