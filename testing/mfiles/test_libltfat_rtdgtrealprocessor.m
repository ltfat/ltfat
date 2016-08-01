function test_libltfat_rtdgtrealprocessor
[~,~,enuminfo]=libltfatprotofile;
LTFAT_FIRWIN = enuminfo.LTFAT_FIRWIN;

gl = 1024;
M = 1500;
a = 200;
W = 10;

gPtr = libpointer('doublePtr',zeros(gl,1));
gdPtr = libpointer('doublePtr',zeros(gl,1));

calllib('libltfat','ltfat_firwin_d',LTFAT_FIRWIN.LTFAT_HAMMING,gl,gPtr);
calllib('libltfat','ltfat_gabdual_painless_d',gPtr,gl,a,M,gdPtr);

plan = libpointer();

calllib('libltfat','ltfat_rtdgtreal_processor_init_win_d',...
    LTFAT_FIRWIN.LTFAT_HAMMING,gl,a,M,W,libpointer('voidPtr'),libpointer('voidPtr'),...
    plan);


[bufIn,fs] = gspi;
bufIn = bsxfun(@times, repmat(bufIn,1,W), [1, rand(1,W-1) + 1]);

%bufIn = ones(1,1000);
bufOut = zeros(size(bufIn));
bufLen = 20;

bufIn = postpad(bufIn,ceil(size(bufIn,1)/bufLen)*bufLen);


for ii=1:length(bufIn)/bufLen
slice = (ii-1)*bufLen + 1 : ii*bufLen;
buf = bufIn(slice,:);
bufInPtr = libpointer('doublePtr',buf);
bufOutPtr = libpointer('doublePtr',zeros(size(buf)));

% Matlab automatically converts Ptr to PtrPtr
calllib('libltfat','ltfat_rtdgtreal_processor_execute_compact_d',plan,bufInPtr,bufLen,W,bufOutPtr);

bufOut(slice,:) = bufOutPtr.Value;

end

inshift = circshift(bufIn,(gl-1));
inshift(1:(gl-1),:) = 0;
plotthat = [bufOut- inshift];
plotthat(end-(gl-1):end) = 0;
stem(plotthat);shg;
%soundsc(bufOut,fs);

calllib('libltfat','ltfat_rtdgtreal_processor_done_d',plan);




 