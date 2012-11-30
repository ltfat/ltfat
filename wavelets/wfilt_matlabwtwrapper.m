function [h,g] = wfilt_matlabwtwrapper(wname)
% wrapper of the Matlab Wavelet Toolbox wfilters function
% remove unnecerary initial zeros of some filters

[lo,hi,lo_s,hi_s] = wfilters(wname);

%     if lo(1)==0 && hi(1)==0 && lo_s(1) == 0 &&  hi_s(1) == 0
%         lo = lo(2:end); hi = hi(2:end);
%         lo_s = lo_s(2:end);hi_s = hi_s(2:end);
%     end

h=cell(2,1);
h{1} = lo;
h{2} = hi;
    

if(nargout>1)
    g=cell(2,1);
    g{1} = lo_s;
    g{2} = hi_s;
end
