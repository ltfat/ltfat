function test_libltfat_rtdgtrealprocessor
[~,~,enuminfo]=libltfatprotofile;
LTFAT_FIRWIN = enuminfo.LTFAT_FIRWIN;

gl = 1024;
M = gl;
a = 100;

gPtr = libpointer('doublePtr',zeros(gl,1));
gdPtr = libpointer('doublePtr',zeros(gl,1));

calllib('libltfat','firwin_d',LTFAT_FIRWIN.HAMMING,gl,gPtr);
calllib('libltfat','gabdual_painless_d',gPtr,gl,a,M,gdPtr);

processorPtr = calllib('libltfat','rtdgtreal_processor_init_d',...
    gPtr,gdPtr,gl,a,M,libpointer('voidPtr'),libpointer('voidPtr'));


[bufIn,fs] = gspi;
%bufIn = ones(1,1000);
bufOut = zeros(size(bufIn));
bufLen = 16;

for ii=1:length(bufIn)/bufLen
slice = (ii-1)*bufLen + 1 : ii*bufLen;
buf = bufIn(slice);
bufInPtr = libpointer('doublePtr',buf);
bufOutPtr = libpointer('doublePtr',zeros(size(buf)));

calllib('libltfat','rtdgtreal_processor_execute_d',processorPtr,bufInPtr,bufLen,bufOutPtr);

bufOut(slice) = bufOutPtr.Value;

end

inshift = circshift(bufIn,(gl-1));
inshift(1:(gl-1)) = 0;
plotthat = [bufOut- inshift];
plotthat(end-(gl-1):end) = 0;
stem(plotthat);shg;
%soundsc(bufOut,fs);

calllib('libltfat','rtdgtreal_processor_done_d',processorPtr);




 