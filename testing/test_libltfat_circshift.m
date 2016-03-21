function test_libltfat_circshift


z = (1:10)+i*(1:10);
zi = complex2interleaved(z);

ziPtr = libpointer('doublePtr',zi);
zoutPtr = libpointer('doublePtr',zi);
calllib('libltfat','circshift_cd',ziPtr,10,2,ziPtr);

ziPtr.Value
interleaved2complex(ziPtr.Value)

