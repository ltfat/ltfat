function test_libltfat_fold

L = 10;
Lfold = 12;
zi = (1:Lfold) ;
zout = zeros(1,Lfold);

ziPtr = libpointer('doublePtr',zi);
zoutPtr = libpointer('doublePtr',zout);


calllib('libltfat','fold_array_d',ziPtr,L,Lfold,zoutPtr);

calllib('libltfat','fold_array_d',ziPtr,L,Lfold,ziPtr);

ziPtr.Value
zoutPtr.Value

%interleaved2complex(zoutPtr.Value)

