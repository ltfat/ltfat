function test_libltfat_fifo
gl = 20;
fifoLen = 120; % Must be at least as big as gl + max expected bufLen
M = gl;
a = 7;
g = firwin('hann',gl);
gd = gabdual(g,a,M);
gg = fftshift(g.*gd)*M;

fifoPtr = calllib('libltfat','rtdgtreal_fifo_init_d',fifoLen,gl,a,1);
fifoPtr.Value.buf.setdatatype('doublePtr',fifoLen+1)
ififoPtr = calllib('libltfat','rtidgtreal_fifo_init_d',fifoLen,gl,a,1);
ififoPtr.Value.buf.setdatatype('doublePtr',fifoLen+gl+1)

bufIn = (1:1000)';
%bufIn = ones(1,1000);
bufOut = zeros(size(bufIn));
bufLen = 100;
bufOutPtr = libpointer('doublePtr',zeros(gl,1));

for ii=1:length(bufIn)/bufLen
slice = (ii-1)*bufLen + 1 : ii*bufLen;
buf = bufIn(slice);
bufInPtr = libpointer('doublePtr',buf);
bufInPtrTmp = libpointer('doublePtr',zeros(size(buf)));

written = calllib('libltfat','rtdgtreal_fifo_write_d',fifoPtr,bufLen,bufInPtr);

while calllib('libltfat','rtdgtreal_fifo_read_d',fifoPtr,bufOutPtr) > 0
    bufOutPtr.Value = bufOutPtr.Value.*gg;
    written = calllib('libltfat','rtidgtreal_fifo_write_d',ififoPtr,bufOutPtr)
end

read = calllib('libltfat','rtidgtreal_fifo_read_d',ififoPtr,bufLen,bufInPtrTmp)

bufOut(slice) = bufInPtrTmp.Value;

end

inshift = circshift(bufIn,(gl-1));
inshift(1:(gl-1)) = 0;
stem([bufOut, inshift]);shg;




calllib('libltfat','rtdgtreal_fifo_done_d',fifoPtr);
calllib('libltfat','rtidgtreal_fifo_done_d',ififoPtr);



 